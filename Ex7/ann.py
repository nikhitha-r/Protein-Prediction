# -*- coding: utf-8 -*-

import os
import sys
import numpy as np

# directory of current script
rel_test_dir = os.path.dirname( __file__ )
# student directory containing the student's solution
rel_root_dir = os.path.split( rel_test_dir )[0] 

if rel_root_dir not in sys.path:
    sys.path.append( rel_root_dir )

class MLP_ANN():
    
    # Possible class labels ( description see below )
    LABELS      = '012HhLU'
    
    # Map all possible labels to an integer (class)
    # During this exercise, we will reduce the problem to a binary decision:
    # Is a given residue part of a transmembrane-helix (1) or not (0)?
    LABEL_2_INT= {  '0' : 0, # non-membrane (inside)
                    '1' : 0, # non-membrane (outside)
                    '2' : 0, # non-membrane (unknown-topology)
                    'L' : 0, # membrane re-entrant loop
                    'U' : 0, # unknown/unresolved
                    'H' : 1, # transmembrane helix
                    'h' : 1  # transmembrane helix
                    }
    
    # get number of classes; this is mainly for experimenting with multi-class 
    NUM_LABELS  = len( set( LABEL_2_INT.values()) ) 
    
    
    ALPHABET     = 'ACDEFGHIKLMNPQRSTVWY' # alphabet of one-letter amino acid codes
    NUM_ALPHABET = len( ALPHABET )        # number of letters in the alphabet


    def __init__( self, nEpochs=501, learning_rate=0.01, half_window_size=3, 
                 n_hidden=50, test_percentage=0.1, seed=42 ):
        '''
            This class is a wrapper for an artifical neural network and all 
            required pre-processing steps. The preprocessing steps include:
                
                1. Store hyperparameters, s.a. number of epochs, learning rate, ...
                2. Retrieve one-hot-encoding for all letters in a given Alphabet
                3. Read in raw data from text file
                4. Encode inputs (one-hot-encoding) and outputs (class labels)
                6. Use 'SAME' padding to artifically add samples at the beginning
                    and the end of the sequence. This ensures that the residues
                    at the start and the end of the sequence can also be predicted
                    by a sliding window approach. Here, zero-padding is used.
                5. Create sliding window view for each residue. As the sequence
                    was already padded in the previous step, the number of residues
                    should not decrease when sliding the window over the sequence.
                6. Split the data set into a train set & a test set
                
            After the preprocessing, the network can be trained. This includes:
                
                - Randomly initialize weight matrices w1, w2. During this step 
                    the bias term has to be taken into account (constant=1)
                - Train for n Epochs by updating your weights after each protein.
                    During this step ONLY the training set is used..
                - The performance of your model is monitored by assessing its
                    performance on your test set. This allows you to check whether
                    your model still learns (train & test loss decrease) or
                    if your model starts to overfit (only train loss decreases).
        '''
        
        self.nEpochs         = nEpochs # number of training epochs
        self.lr              = learning_rate # the learning rate
        # window size, e.g. 7= xxxXxxx, with X being the predicted residue and 
        # x being neighbouring residues; giving half window size allows only uneven numbers
        self.half_window     = half_window_size # half window size
        self.window_size     = self.half_window * 2 + 1 # window size
        self.H               = n_hidden # number of nodes in the hidden layer
        self.test_percentage = test_percentage # percentage of samples used for testing
        self.seed            = seed
        self.bias            = 1 # bias term. 
        
        file_dir = os.path.dirname( __file__ ) # path to data set (text file)
        self.path            = os.path.join( file_dir, *['tests', 'opm_tmps.fasta'] )
        
        # input size for each residue = 
        # number of unique amino acid letters (20) * window size + bias (1)
        D_in  = self.NUM_ALPHABET * self.window_size + self.bias 
        D_out = self.NUM_LABELS # output dimension; number of predicted labels
        
        # Randomly initialize weights for input->hidden (w1) and hidden->output (w2) layer
        rnd_generator = self.__reinit_rnd_state() 
        self.w1  = rnd_generator.randn( D_in,               self.H )
        self.w2  = rnd_generator.randn( self.H + self.bias, D_out  )
    
    ######################## UTILITY FUNCTIONS ################################
    
    def __reinit_rnd_state( self ):
        '''
            Re-initialize the random state generator every time you need random
            variables in different functions. This avoids that your tests fail 
            just because you are completing the exercise sheet in a different 
            order as we expect it.
        '''
        return np.random.RandomState( self.seed )
    
    def get_window_size( self ):
        # return window size. Required for testing. Do not change
        return self.window_size
    
    def get_alphabet_size( self ):
        # return alphabet size. Required for testing. Do not change
        return self.NUM_ALPHABET
    
    
    ########################### PREPROCESSING #################################
    
    def _get_raw_data( self ):
        '''
            This function reads in the raw data from a text-file given by
            self.path. The file format should already be known from exercise 6:
                > HEADER
                    SEQUENCE
                    TRANSMEMBRANE-ANNOTATION
            Create two dictionaries, one for inputs, another one for outputs.
            The keys for these dictionaries are the unique identifiers given in
            the header; e.g. P61829. The identifier is always the first identifier
            after the '>' and before the first '|'. These identifiers are used
            throughout this exercise.
            The values are strings coding for the sequence (inputs) or the 
            corresponding transmembrane annotation (targets).
            Return the two dictionaries as a tuple.
            
            Inputs:
                    self
            Returns:
                    Tuple of two dictionaries, summarizing the inputs and the
                    outputs (in this order) with identifiers as keys and strings
                    as values.
        '''
        
        inputs = dict() # create containers for inputs
        outputs= dict() # and outputs
        
        with open( self.path, 'r' ) as data_file: # open file handle
            for line in data_file: # for each line in the file
                # TODO: Read in raw data
                #print(line)
                if ">" in line:
                    #start = s.index( first ) + len( first )
                    end = line.index( "|", 0 )
                    key = line[1:end]
                elif key not in inputs:
                    inputs[key] = line.replace("\n", "")
                else:
                    outputs[key] = line.replace("\n", "")
        self.raw_inputs  = inputs
        self.raw_outputs = outputs
        return (inputs, outputs)
    
    
    def _get_one_hot_encoding( self ):
        '''
            This function encodes all letters given in ALPHABET (see above) as
            one-hot encoded vectors. Here, it is important to maintain the 
            order of the original ALPHABET in order to use the same encoding
            as we do during testing your results.
            E.g. A -> [ 1, 0, 0, ..]; C -> [ 0, 1, 0, ..]; D -> [ 0, 0, 1, ..]
        
            Inputs:  self (ALPHABET is accessible via self.)
            Returns: Dictionary which holds amino acid letters as keys and the
                     corresponding one-hot-encoding as numpy array.
        '''
        one_hot_dict     =  {}
        # TODO: Implement one-hot-encoding dictionary
        zero_list = [0 for j in range(self.NUM_ALPHABET)]
        for i in range(self.NUM_ALPHABET):
            
            if self.ALPHABET[i] not in one_hot_dict:
                new_list = zero_list.copy()
                new_list[i] = 1
                one_hot_dict[self.ALPHABET[i]] = new_list
        self.AA_2_INT    = one_hot_dict
        return one_hot_dict
    
    
    def _get_one_hot_inputs( self ):
        '''
            Before training the network you have to transform your input data to 
            make it usable/understandable for the neural network.
            For this we are using the one-hot-encoding created by your function
            '_get_one_hot_encoding' which is accessible via 'self.AA_2_INT'.
            Iterate over all sequence in your input data and transform every
            residue to the corresponding one-hot-vector.
            Store the result again as dictionary, with the original identifiers
            being again the keys and the one-hot-encoded sequences being the 
            values ( numpy arrays ).
            The one-hot-encoded sequence have to have this shape: 
            ( n_rows, n_columns ) = ( length_of_protein, length_of_Alphabet )
            
            Inputs:
                self ( the raw data should already be read in by your function 
                      '_get_raw_data' by now and stored in self.raw_inputs)
            Returns:
                dictionary, holding one-hot-encoded sequences as numpy arrays 
                (values) and identifiers as keys.
        '''
        
        encoded_inputs=dict()
        
        for uniprot_id, nominal_data in self.raw_inputs.items(): # for each sample
            # TODO: Transform inputs via one-hot-encoding [self.AA_2_INT]
            encod_val = np.zeros((len(nominal_data), self.NUM_ALPHABET))
            for i in range(np.shape(encod_val)[0]):
                encod_val[i][:] = self.AA_2_INT[nominal_data[i]]
            encoded_inputs[uniprot_id] = encod_val
            
        self.encoded_inputs = encoded_inputs
        print("######################")
        return encoded_inputs


    def _get_sliding_window_view( self ):
        '''
            This function iterates over all training samples (one-hot-encoded)
            in order to create a sliding window view for each residue. Every
            residue will be represented by its own one-hot-encoding as well as
            the one-hot-encoding of all neighbors within window_size.
            While creating the sliding window view, the ordering of residues does
            not change; meaning, that the actual residue of interest is in the
            middle of your sliding window while the neighboring residues are added
            to the left or right.
            The size of the sliding window is defined during class initialization.
            Use zero-vectors for padding the sequences in such a way that the 
            length of your sequence does not decrease when sliding the window
            over the sequence.
            
            Inputs:
                self (encoded inputs should be available via self.encoded_inputs)
            Returns:
                dictionary with protein identifiers as keys and numpy array as
                values. The shape of the numpy array is now:
                ( n_rows, n_columns ) = ( length_of_protein, window_size * len(ALPHABET))
        '''
        
        sliding_window_inputs = dict()
        padding_size          = self.half_window
        print(self.window_size)
        for uniprot_id, encoded_input in self.encoded_inputs.items(): # for each sample
            # TODO: Implement sliding window view for each residue
            final_padded_input = np.zeros((encoded_input.shape[0], encoded_input.shape[1] * self.window_size))
            padded_input = np.zeros((encoded_input.shape[0] + 2*padding_size, encoded_input.shape[1]))
            padded_input[padding_size:padding_size+encoded_input.shape[0], :] = encoded_input
            for i in range(padding_size, encoded_input.shape[0]+ padding_size):
                old_encod = padded_input[i, :]
                #new_encod = np.array((encoded_input.shape[1] * self.window_size))
                new_encod = np.array([])
                for j in range(i-padding_size, i):
                    new_encod = np.hstack([new_encod, padded_input[j]])
                new_encod = np.hstack([new_encod, padded_input[i]])
                for j in range(i+1, i+padding_size+1):
                    new_encod = np.hstack([new_encod, padded_input[j]])
                final_padded_input[i-padding_size,:] = new_encod
                #print(new_encod)
            print(final_padded_input.shape)
            np.savez("f.npz", final_padded_input)
            sliding_window_inputs[uniprot_id] = final_padded_input
            #print(final_padded_input)
            
            
            
        self.inputs = sliding_window_inputs
        return sliding_window_inputs
    
    
    def _get_integer_outputs( self ):
        '''
            This function maps all target labels to their corresponding classes.
            Use the static variable LABEL_2_INT (see above) for this mapping.
            Iterate over all items in self.raw_outputs (dictionary created in 
            get_raw_outputs), transform String labels to integer classes and 
            save these as numpy array.
            While doing so, fill a second dictionary of the same structure 
            (keys are identifiers, values are numpy arrays of the same size as
            the protein sequence) which holds information on whether a given 
            residue in a protein is resolved (1) or not (0).
            The intuition behind this is that those residues can be masked out 
            during training, meaning that we remove those residues as we do not 
            know whether they are part of a transmembrane helix or not.
            
            Inputs:
                self
            Returns:
                Tuple of two dictionaries, one holding data on targets, the other 
                holding data on which residues get masked out during training.
                Keys: Identifiers, Values: Numpy arrays
        '''
        
        encoded_outputs = dict()
        mask_dict        = dict()

        for uniprot_id, nominal_data in self.raw_outputs.items(): # for each sample
            # TOOD: Transform class labels to integers
            label_data = np.array([])
            mask_data = np.array([])
            for i in range(len(nominal_data)):
                label_data = np.hstack([label_data, self.LABEL_2_INT[nominal_data[i]]])
                if nominal_data[i] == "U":
                    mask_data = np.hstack([mask_data, 0])
                else:
                    mask_data = np.hstack([mask_data, 1])
            encoded_outputs[uniprot_id] = label_data
            mask_dict[uniprot_id] = mask_data
            
        self.outputs = encoded_outputs
        self.mask    = mask_dict
        return ( encoded_outputs, mask_dict )
    
    
    def _remove_unresolved( self ):
        '''
            Iterate over the sequences in your data set and remove those residues
            which have the label 'U' (unresolved). These residues were not resolved
            in the structure so we can't say whether they are part of the 
            transmembrane helix or not. However, we still know which type of 
            amino acid is at this position, so were able to use these residues
            in the sliding window of neighboring residues. Therefore, the order
            of pre-processing steps is important.
            Use the self.mask data structure which you've created in 
            '_get_integer_outputs' for removing unresolved residues from
            self.inputs & self.outputs.
            Inputs:
                self ( inputs, outputs and mask should be available by now)
            Returns:
                Inputs & outputs without unresolved residues
        '''

        for uniprot_id, input_data in self.inputs.items():
            # TODO: Remove unresolved residues
            mask_data = self.mask[uniprot_id]
            rem_indxs = [i for i in range(len(mask_data)) if mask_data[i] == 0]
            if len(rem_indxs) > 0:
                input_data = np.delete(input_data, rem_indxs, axis=0)
                self.inputs[uniprot_id] = input_data
                output_data = self.outputs[uniprot_id]
                output_data = np.delete(output_data, rem_indxs, axis=0)
                self.outputs[uniprot_id] = output_data
                
            
            
        return ( self.inputs, self.outputs )
    
    
    def _get_train_test_split( self ):
        '''
            This function splits the proteins in your one-hot-encoded inputs 
            and your integer-encoded outputs into a training and a test set.
            The percentage size of your test set is defined by the 
            'test_percentage' variable defined during object initialization.
            You have to return one tuple which contains four dictionaries:
                train_inputs & train_outputs hold proteins used for training
                test_inputs  & test_outputs  hold proteins used for testing
            Split the data set only on the level of proteins, not residues.
            
            Inputs:
                self
            Returns:
                Tuple containing four dictionaries, 
                two for training (inputs & outputs) and two for testing.
                The data structure of these dictionieres does not change;
                keys are always identifiers, and values are numpy arrays either
                containing one-hot-encoded inputs or integer-encoded classes.
        '''
        train_inputs = dict()
        train_outputs= dict()
        test_inputs  = dict()
        test_outputs = dict()

        # init a new generator object for random numbers
        rnd_generator = self.__reinit_rnd_state() 
        
        for uniprot_id, input_data in self.inputs.items(): # for every sample
            
            rnd         = rnd_generator.rand()
            output_data = self.outputs[ uniprot_id ] # get output
            
            if rnd < self.test_percentage: # add samples to test set
                # TODO: Add sample to test dictionaries (input & output)
                test_inputs[uniprot_id] = self.inputs[uniprot_id]
                test_outputs[uniprot_id] = self.outputs[uniprot_id]
            else: # add sample to train set
                
                # TODO: Add sample to train dictionaries (input & output)
                train_inputs[uniprot_id] = self.inputs[uniprot_id]
                train_outputs[uniprot_id] = self.outputs[uniprot_id]

        self.train_inputs  = train_inputs
        self.train_outputs = train_outputs 
        self.test_inputs   = test_inputs
        self.test_outputs  = test_outputs
        return ( train_inputs, train_outputs, test_inputs, test_outputs )
    


    def get_inputs_and_outputs( self ):
        '''
            Wrapper function which summarizes all functions required for the
            pre-processing of inputs & outputs.
            Call this function before training the network in order to make
            inputs & outputs accessible.
            This function does not have to return anything as the pre-processed
            inputs & outputs are stored in class variables self.train_inputs,
            self.train_outputs, self.test_inputs and self.test_outputs.
        '''

        # read in raw data from text file and store inputs and outputs in dict
        self._get_raw_data()
        
        # one hot encoded single letter amino acid codes
        self._get_one_hot_encoding() 
        
        # Encode input sequences via one-hot-encoding 
        self._get_one_hot_inputs()
        
        # Create sliding window view for each residue
        self._get_sliding_window_view()
        
        # Encode class labels as integers and save unresolved residues as mask
        self._get_integer_outputs()
        
        # remove residues which were not resolved in the structure from the set.
        # Please understand that this has to be done after creating the sliding
        # window view as these residues appear in the sequence but could not be
        # resolved during solving the structure.
        self._remove_unresolved()
        
        # Split proteins in a training and a test set
        self._get_train_test_split()
        
        return None


    ########################## DATA ANALYSIS #################################
    def get_num_samples( self ):
        '''
            This function returns the number of all samples ( single residues 
            without unresolved residues ) in your data set.
            Inputs:
                self ( self.inputs or self.outputs should be available by now)
            Returns:
                Number of samples (residues) in the whole data set
        '''
        n_samples = 0.
        # TODO: Return the number of samples
        for uniprot_id, input_data in self.inputs.items():        
            n_samples += input_data.shape[0]    
        return n_samples
    
    
    def get_num_masked_out( self ):
        '''
            This function returns the number of samples (single residues) in 
            your data set which were unresolved (label='U').
            Inputs:
                self ( self.mask should be available by now )
            Returns:
                Number of samples (residues) which are masked out during training
        '''
        n_masked_out = 0.
        # TOOD: Count unresolved samples
        for uniprot_id, input_data in self.inputs.items():
            # TODO: Remove unresolved residues
            mask_data = self.mask[uniprot_id]
            rem_indxs = [i for i in range(len(mask_data)) if mask_data[i] == 0]
            n_masked_out += len(rem_indxs)
        return n_masked_out

    
    def get_num_pos_and_neg( self ):
        '''
            This function returns the number of positive samples ( residues in 
            a transmembrane helix) and negative samples (not part of TMH) as a
            tuple. Calculate the values for the modified outputs, e.g. not 
            containing unresolved residues.
            Inputs:
                self ( self.outputs should be available by now )
            Returns:
                Tuple consisting of (num_pos, num_neg)
        '''
        n_pos = 0.
        n_neg = 0.
        # TODO Count positive (transmembrane) and negative samples
        for uniprot_id, output_data in self.outputs.items():
            for i in range(len(output_data)):
                if output_data[i] == 0:
                     n_neg += 1
                else:
                    n_pos += 1
            
        return ( n_pos, n_neg )
        
    
    ############################## TRAINING  ##################################    
    
    def _stable_softmax( self, predictions ):
        """
        Calculate the softmax 'normalization' for raw prediction scores.
        'Normalization', because the raw prediction scores are transformed to
        a probability within [0,1]. Subtract a constant value of max(predictions)
        to the exponents to make the softmax calculation numerically stable.
        This constant factor does not affect our training progress as we 
        multiply with a constant factor (learning rate) anyway.
        inputs:
            self
            predictions: numpy array of predicted lables for a whole batch 
                            (here: protein)
        Output:
            Normalized predictions for one batch (protein) as numpy array
            Expected shape: length of protein x 2
        """
        softmax = np.zeros((predictions.shape))
        max_pred = np.max(predictions)
        
        for i in range(predictions.shape[0]):
            vals = []
            for j in range(predictions.shape[1]):
                vals.append(np.exp(predictions[i, j] - max_pred))
            vals = vals/np.sum(vals)
            softmax[i, :] = vals
        # TODO: Implement softmax 
        return softmax
    
    
    def _cross_entropy( self, y_pred, y ):
        """
            Calculate the cross entropy loss for a given prediction (y_pred) and
            the corresponding true class label (y).
            Please understand that in the context of machine learning (when
            calculating error rates between 0 and 1), the cross_entropy loss
            is equal to the log loss. Also double check the information on the
            slides if necessary.
            Hint: When calculating the loss, think of the groundtruth labels in
                    y as an array of indices. Also think about a
                    possible simplification of the cross-entropy summation if 
                    one of the labels is always encoded via 0 (multiplication
                    with 0).
            Inputs:
                self
                y_pred: the probability distribution of predicted class labels
                        Shape: ( num_examples x num_classes )
                y     : the true class label [0,1].
                        Shape: ( num_examples x 1 )
            Output:
                loss: Difference between predicted and true class label.
                Again this is calculated for a whole batch (protein). However,
                the loss for all samples in a single batch is simply summed at 
                the end. The expected return value is hence a single scalar.
        """
        loss = 0.
        # TODO: Implement cross entropy loss
   
        for i in range(y.shape[0]):
            if y[i] != 0:
                calc_loss = y[i] * np.log(y_pred[i][1])
                loss += np.sum(calc_loss)
            else:
                calc_loss = np.log(y_pred[i][0])
                loss += calc_loss
        loss = -loss
      

        return loss
    
    
    def _delta_cross_entropy( self, y_pred, y ):
        """
            Calculate the derivative of the cross entropy loss (or log loss) 
            with respect to y_pred. 
            This gives you the gradient for changing your weigths.
            Please understand that in the context of machine learning (when
            calculating error rates between 0 and 1), the cross_entropy loss
            is equal to the log loss. Also double check the information on the
            slides if necessary.
            
            Hint 1:
                Consider the batch size when returning the gradient of the loss by
                dividing the gradients by the number of samples in the batch.
            Hint 2:
                While calculating the derivative of cross entropy, 
                consider that you have to take the derivative of the softmax into 
                account as well. This allows to simplify the final expression.
            Hint3: 
                When calculating the delta, think of the groundtruth labels in
                y rather as a matrix with a one-hot encoding of the label 
                ([1, 0] or [0, 1]) than of an array. Also think about a
                possible simplification of the cross-entropy summation if 
                one of the labels is always encoded via 0 (multiplication w. 0).

            Inputs:
                self
                y_pred: the probability distribution of predicted class labels
                        Shape: ( num_examples x num_classes )
                y     : the true class label.
                        Shape: ( num_examples x 1 )
            Output:
                gradient: The gradient for changing the weights. Again, remember
                that we are performing batch learning, meaning that we do not 
                update the weights after each sample separately but that we 
                accumulate the loss of multiple samples and update the weights 
                simultaneously after all samples in one batch are processed.
                Here, you should return a numpy array of shape 
                    ( length of the protein x 2)
                You can use the dot product to update weights simultaneously
                during backpropagation afterwards.
        """
        gradient = np.zeros((y_pred.shape))
        # TODO: Implement gradient of cross-entropy loss
       
        num_sam = len(y)
        for i in range(len(y)):
            if y[i] == 1:
                y_new = [0, 1]
            else:
                y_new = [1, 0]
           
            #print(y_pred[i])
            gradient[i, :] = (y_pred[i] - y_new)/num_sam
        return gradient
    

    def train_network( self ):
        '''
            Wrapper function for training (iterating over all epochs & samples)
            and testing.
        '''

        for epoch in range( self.nEpochs ): # for each epoch
            # for each sample in the training sample
            for uniprot_id, x_train in self.train_inputs.items(): 
                
                y_train      = self.train_outputs[ uniprot_id ] # get outputs
                self._predict( x_train, y_train )         # train 
            
            if epoch%20==0: # monitor loss
                acc_train, avg_loss_train, conf_mat_train = self._test_performance( 
                                                                eval_train_set = True  )
                acc_test,  avg_loss_test,  conf_mat_test  = self._test_performance( 
                                                                eval_train_set = False )
            
        return self.w2
        
    
    def _predict( self, x, y, is_training=True ):
        '''
            Function for predicting one batch (here: one protein) of samples.
            If the boolean 'is_training' is set to True, the weights are adjusted
            according to the loss. The predictions are always returned.
            Using a boolean here allows us to use one function for training
            and testing.
            Again, please remember that we are performing batch learning while
            training the network. This means that we do not update the weights
            after each single samples but after all samples in one batch
            simultaneously. The gradient which is returned should contain the
            gradient for every single sample in a batch (len_protein x 2) and 
            by using the dot product you can update all weights simultaneously 
            for all samples.
            Also remember to add a bias term to the input layer as well as the
            hidden layer.
            Inputs:
                self
                x: One batch (protein) of one-hot-encoded Inputs (sliding window)
                y: One batch (protein) of integer-encoded labels
                is_training: Boolean, controlling whether backpropagation (True)
                    is used or not (False).
            Output:
                y_pred: Predicted labels for this batch
        '''
        # Forward pass: compute predicted y
        bias = np.ones((x.shape[0], 1)).tolist()
        x = np.column_stack((x, bias))
        
        h1 = np.dot(self.w1.T, x.T)
        #print(h1)
        h_relu = np.maximum(0, h1).T
        #print("ad", h1)
        
        bias = np.ones((h_relu.shape[0], 1)).tolist()
        h_relu = np.column_stack((h_relu, bias))
       
        y_out = np.dot(h_relu, self.w2)
        y_pred = self._stable_softmax(y_out)
        # Backprop to compute gradients of w1 and w2 with respect to loss
        if is_training: # during training, change the weights based on loss 
            grad_y_pred   = self._delta_cross_entropy( y_pred, y  )
            """
            print(grad_y_pred)
            print(grad_y_pred.shape)
            print(h_relu.shape)
            """
            # Update weights only while training, not testing
            dw2 = np.dot(h_relu.T, grad_y_pred) 
            """
            print("***", dw2.shape)
            print(self.w2.shape)
            print(self.w1.shape)
            """
            grad_hidden = np.dot(grad_y_pred, self.w2[0:-1].T)
            #print("**********")
            #print(grad_hidden.shape)
            grad_h = grad_hidden.T.copy()
            #print(np.where(h1 < 0))
            grad_h[np.where(h1 < 0)] = 0
           
            #print(x.shape)
            #print(grad_h.shape)
            dw1 = np.dot(grad_h, x).T
            """
            print(dw1.shape)
            print(self.w1.shape)
            """
            self.w1 -= self.lr * dw1
            self.w2 -= self.lr * dw2
            
            
        return y_pred
    
    
    def _test_performance( self, eval_train_set ):
        '''
            This function is used to monitor the performance throughout the 
            training progress. The difference in performance on the training set 
            and the test set allows you to stop training before you overfit your 
            model.
            
            Input:
                self
                eval_train_set: Boolean, describing whether the current model
                    should be evaluated using the training set or the test set.
            Output:
                Tuple summarizing several performance values used for monitoring
                the training progress.
                The tuple should consists of ( accuracy, average loss and 
                confusion matrix ).
                Accuracy         = Number of correctly predicted samples / 
                                   Number of all samples in this set
                Average loss     = Average of loss over all batches
                Confusion matrix = Rows are true labels and columns are predicted
                                    labels:
                                    TN = True  Negatives | FP = False Positives
                                    FN = False Negatives | TP = True  Positives
        '''
        
        if eval_train_set:
            inputs  = self.train_inputs
            outputs = self.train_outputs
        else:
            inputs  = self.test_inputs
            outputs = self.test_outputs
        
        # reset statistics
        conf_mat = np.zeros( (self.NUM_LABELS, self.NUM_LABELS), dtype=np.int )
        avg_loss = 0. # average loss
        total    = 0. # total number of residues
        correct  = 0. # number of correctly predicted samples
        acc      = 0. # prediction accuracy
        fp = 0
        fn = 0
        tn = 0
        tp = 0
        loss = 0
        # TODO: Implement performance estimates
        #print(inputs)
        for uniprot_id, x_train in inputs.items(): 
            y_pred = self._predict(inputs[ uniprot_id ], outputs[ uniprot_id ], is_training=False)
            gt = outputs[uniprot_id]
            loss += self._cross_entropy( y_pred, gt )
            for i in range(len(gt)):
                pred = np.argmax(y_pred[i])
                if gt[i] == 0 and pred == 0:
                    tn += 1
                elif gt[i] == 0 and pred == 1:
                    fp += 1
                elif gt[i] == 1 and pred == 0:
                    fn += 1
                elif gt[i] == 1 and pred == 1:
                    tp += 1
        n_samples = 0
        correct = tn + tp
        for uniprot_id, input_data in inputs.items():        
            n_samples += input_data.shape[0]  
        acc = correct/n_samples
        avg_loss = loss/len(inputs)
        conf_mat = np.array([[tn , fp], [fn, tp]])
        print(conf_mat)

        return acc, avg_loss, conf_mat
    
def main():
    print("d")
    
    ann_numpy = MLP_ANN()
    ann_numpy.get_inputs_and_outputs()
    ann_numpy.train_network()
    

if __name__ == '__main__':
    main()
"""
import os
import json
relative_path = "/Users/nikhitha/Documents/Protein Prediction/Exercises/pp1cs2020exercise7-ge73tag/tests"
with open(os.path.join(relative_path, 'test.json')) as json_file:
    json_data = json.load(json_file)
bias = np.load(os.path.join(relative_path, json_data["sliding_window_inputs"]))
s = bias
out = np.load(os.path.join("/Users/nikhitha/Documents/Protein Prediction/Exercises/pp1cs2020exercise7-ge73tag/f.npz"))
o =out
"""