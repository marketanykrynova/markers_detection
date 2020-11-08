function [entropy_of_words_for_given_word_lengths, sequence_without_ambiguous_bases_and_modified_gaps] = calculate_entropy(length_of_reference, length_of_words)
% [ENTROPY_OF_WORDS_FOR_GIVEN_WORD_LENGTHS, MODIFIED_SEQUENCES] = CALCULATE_ENTROPY(LENGTH_OF_REFERENCE, LENGTH_OF_WORDS) Run the program "Calculate entropy of words"
%   length_of_reference: length of the sequence used for reference-based assembly
%   length_of_words: length of words used for entropy calculation
%   entropy_of_words_for_given_word_lengths: structure that contains value of entropy of words for given lengths
%   sequence_without_ambiguous_bases_and_modified_gaps: sequences of isolates after preprocessing
%   example: [entropy_of_words_for_given_word_lengths, sequence_without_ambiguous_bases_and_modified_gaps] = calculate_entropy(2955294, [50 100 150 200 250 300 350 400])

[headers, sequences] = load_sequences(length_of_reference);

save('headers.mat','headers')

sequence_without_ambiguous_bases = replace_ambiguous_nucleotides(sequences);
sequence_without_ambiguous_bases_and_modified_gaps = modified_gaps(sequence_without_ambiguous_bases);
entropy_of_words_for_given_word_lengths = words_entropy(length_of_words,sequence_without_ambiguous_bases_and_modified_gaps);

save('sequence_modified_ambiguous_nt_and_gaps.mat','sequence_without_ambiguous_bases_and_modified_gaps')

function [headers, sequences] = load_sequences(length_of_reference)
% [HEADERS, SEQUENCES] = LOAD_SEQUENCES(LENGTH_OF_REFERENCE) Load all FASTA sequences in the folder.
%   length_of_reference: length of the sequence used for reference-based assembly
%   headers: headers of sequences from FASTA files
%   sequences: aligned nucleotide sequences from FASTA files

list = dir('*.fasta');

if ~isempty(list)
    for i = 1:length(list)
        
        [headers{i},seq{i}] = fastaread(list(i).name);
        
        if length(seq{i}) ~= length_of_reference
            diff = length_of_reference - length(seq{i});
            seq{i} = [seq{i} repmat('-',1,diff)];
        end
        
        sequences(i,:) = seq{i};
        
    end
else
    error('There are no FASTA files in current folder.')
end

function sequence_without_ambiguous_bases = replace_ambiguous_nucleotides(sequences)
% SEQUENCES_WITHOUT_AMBIGUOUS_BASES = REPLACE_AMBIGUOUS_NUCLEOTIDES(SEQUENCES) Preprocessing of aligned sequences - removing ambiguous nucleotides and replacing them by gaps
%   sequences: aligned nucleotide sequences from FASTA files
%   sequences_without_ambiguous_bases: aligned sequences without ambiguous nucleotides

amb_positions = sum(sequences == 'R' | sequences == 'Y' | sequences == 'S' | sequences == 'W' | sequences == 'K' | sequences == 'M');
sequence_without_ambiguous_bases = sequences;
sequence_without_ambiguous_bases(:,amb_positions>0) = '-';

function sequence_without_ambiguous_bases_and_modified_gaps = modified_gaps(sequence_without_ambiguous_bases)
% SEQUENCES_WITHOUT_AMBIGUOUS_BASES_AND_MODIFIED_GAPS = REPLACE_AMBIGUOUS_NUCLEOTIDES(SEQUENCES_WITHOUT_AMBIGUOUS_BASES) Preprocessing of aligned sequences - removing ambiguous nucleotides and replacing them by gaps
%   sequences_without_ambiguous_bases: aligned sequences without ambiguous nucleotides
%   sequences_without_ambiguous_bases_and_modified_gaps: aligned sequences without ambiguous nucleotides and with modified gaps

gaps = sum(sequence_without_ambiguous_bases == '-');
sequence_without_ambiguous_bases_and_modified_gaps = sequence_without_ambiguous_bases;
sequence_without_ambiguous_bases_and_modified_gaps(:,gaps>0) = '-';

function entropy_of_words_for_given_word_lengths = words_entropy(length_of_words, sequence_without_ambiguous_bases_and_modified_gaps)
% [ENTROPY_OF_WORDS_FOR_GIVEN_WORD_LENGTHS] = WORDS_ENTROPY(LENGTH_OF_WORDS, SEQUENCES_WITHOUT_AMBIGUOUS_BASES_AND_MODIFIED_GAPS) Calculate the values of words entropy and save it to structure and to separated mat files
%   length_of_words: length of words used for entropy calculation
%   sequences_without_ambiguous_bases_and_modified_gaps: aligned sequences without ambiguous nucleotides and with modified gaps
%   entropy_of_words_for_given_words_length: structure that contains value of entropy of words for given lengths

entropy_of_words_for_given_word_lengths = struct();

number_of_seq = size(sequence_without_ambiguous_bases_and_modified_gaps,1);
length_of_seq = size(sequence_without_ambiguous_bases_and_modified_gaps,2);

for num_of_window = 1:length(length_of_words)
    
    a = 1;
    
    for i = 1:length_of_seq - length_of_words(num_of_window) + 1
        
        seq = sequence_without_ambiguous_bases_and_modified_gaps(:,i:i + length_of_words(num_of_window) - 1);
        [unique_words, ~, cluster] = unique(seq,'rows');
        
        for j = 1:size(unique_words,1)
            number_of_words_in_each_cluster(j) = sum(cluster == j);
        end
        
        
        entropy_of_words(a) = -1*(sum((number_of_words_in_each_cluster./number_of_seq).*(log2(number_of_words_in_each_cluster./number_of_seq))));
        
        a = a+1;
        clearvars  number_of_words_in_each_cluster
        
    end
    
    name_of_field = ['word_length_' num2str(length_of_words(num_of_window))];
    entropy_of_words_for_given_word_lengths.(name_of_field) = entropy_of_words;
   
    filename = ['entropy_of_words_window_' num2str(length_of_words(num_of_window)) '.mat'];
    save(filename,'entropy_of_words')
    clearvars seq entropy_of_words
    
end



