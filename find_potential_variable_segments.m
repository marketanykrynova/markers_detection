function results = find_potential_variable_segments(genbank_file, threshold_level)
% RESULTS = FIND_POTENTIAL_VARIABLE_SEGMENTS(GENBANK_FILE, THRESHOLD_LEVEL) Run the program "Find potential variable segments"
%   genbank_file: file name of genbank file with sequence used for reference-based assembly
%   threshold_level: value that is used to set threshold
%   results: structure with variable markers
%   example: results = find_potential_variable_segments('NC_017022.1.gb', 0.4)

% load files created during 'calculate_entropy.m'
load ('sequence_modified_ambiguous_nt_and_gaps.mat');
load ('headers.mat');
list = dir('entropy_of_words_window_*');

number_of_seqs = size(sequence_without_ambiguous_bases_and_modified_gaps,1);
results = struct();

for num_of_window = 1:length(list)
    
    num = regexp(list(num_of_window).name, '\d+', 'match');
    window = str2double(num{1});
    
    load (list(num_of_window).name);
    
    peaks = get_coordinates_of_potential_var_parts(window, threshold_level, entropy_of_words);
    [all_seq, field_names] = find_one_variable_segment(peaks, sequence_without_ambiguous_bases_and_modified_gaps, headers, number_of_seqs);
    [mat, num_of_clusters_two_seqs] = find_combination_of_two_variable_segments(all_seq, headers, number_of_seqs, field_names);
    results = save_potential_var_parts_to_results(results, mat, num_of_clusters_two_seqs, window, peaks, all_seq, field_names);
    
    clearvars peaks mat all_seq

    
end

results(1) = [];
results = get_information_about_segments_from_gb_file(genbank_file, results);
results = get_only_results_with_max_number_of_clusters(results);

save('results.mat', 'results')

function peaks = get_coordinates_of_potential_var_parts(window, threshold_level, entropy_of_words)
% PEAKS = GET_COORDINATES_OF_POTENTIAL_VAR_PARTS(WINDOW, THRESHOLD_LEVEL, ENTROPY_OF_WORDS) x
%   window:
%   threshold_level:
%   entropy_of_words:
%   peaks:

threshold = max(entropy_of_words) - threshold_level;
[~,y] =  find(entropy_of_words>=threshold);

if isempty(y)
    return
end

difference = diff(y);
peaks(1,1) = y(1);

i = 1;
row = 1;

if length(y) == 1
    peaks(1,2) = y(1) + window - 1;
end

while i <= length(difference)
    
    if difference(i) <= window
        peaks(row,2) = y(i+1) + window - 1;
        i = i + 1;
    else
        peaks(row,2) = y(i) + window - 1;
        i = i + 1;
        row = row + 1;
        peaks(row,1) = y(i);
    end
    
end

if peaks(row,2) == 0
    peaks(row,2) = peaks(row,1) + window - 1;
end

function [all_seq, field_names] = find_one_variable_segment(peaks, sequence_without_ambiguous_bases_and_modified_gaps, headers, number_of_seqs)
% [ALL_SEQ, FIELD_NAMES] = FIND_ONE_VARIABLE_SEGMENT(PEAKS, SEQUENCE_WITHOUT_AMBIGUOUS_BASES_AND_MODIFIED_GAPS, HEADERS, NUMBER_OF_SEQS) Analyze potential variable fragment and determine number of clusters from phylogenetic tree.
%   peaks: peaks of entropy signals which are above the threshold
%   sequence_without_ambiguous_bases_and_modified_gaps: aligned sequences without ambiguous nucleotides and with modified gaps
%   headers: headers of sequences from FASTA files
%   number_of_seqs: number of isolates
%   all_seq: sequences which corresponded to peaks
%   field_names: names of sequences which corresponded to peaks

num_of_clusters = [];

for i = 1:size(peaks,1)
        seq = sequence_without_ambiguous_bases_and_modified_gaps(:,peaks(i,1):peaks(i,2));
        d = seqpdist(seq, 'Alphabet', 'NT','Indels','complete-del', 'Method','Kimura');
    try
        tree = seqlinkage(d, 'average', headers);
        clusters = cluster(tree,[],'criterion','maximum','maxclust',number_of_seqs);
        num_of_clusters(end+1) = length(unique(clusters));
        name = ['seq_' num2str(i)];
        all_seq.(name) = seq;
    catch
        disp('Tree cannot be constructed.')
        num_of_clusters(end+1) = 0;
    end
    
end

field_names = fieldnames(all_seq);

if length(field_names) == 1
    continue
end

function [mat, num_of_clusters_two_seqs] = find_combination_of_two_variable_segments(all_seq, headers, number_of_seqs, field_names)
% [MAT, NUM_OF_CLUSTERS_TWO_SEQS] =
% FIND_COMBINATION_OF_TWO_VARIABLE_SEGMENTS(ALL_SEQ, HEADERS, NUMBER_OF_SEQS, FIELD_NAMES) Find combination of two variable fragments and number of clusters.
%   all_seq: sequences which corresponded to peaks
%   headers: headers of sequences from FASTA files
%   number_of_seqs: number of isolates
%   field_names: names of sequences which corresponded to peaks
%   mat: matrix with number of clusters for combination of two variable segments
%   num_of_clusters_two_seqs: vector with number of clusters for combination of two variable segments

num_of_clusters_two_seqs = [];

for i = 1:length(fieldnames(all_seq)) - 1
    for j = 2:length(fieldnames(all_seq))
        if j ~= i && i<j
            
            two_seqs = [all_seq.(field_names{i}) all_seq.(field_names{j})];
            try
                d = seqpdist(two_seqs, 'Alphabet', 'NT','Indels','complete-del', 'Method','Kimura');
                tree = seqlinkage(d, 'average', headers);
                clusters = cluster(tree,[],'criterion','maximum','maxclust',number_of_seqs);
                num_of_clusters_two_seqs(end+1) = length(unique(clusters));
                mat(i,j) = length(unique(clusters));
            catch
                disp('Tree cannot be constructed.')
                mat(i,j) = 0;
            end
        end
    end
end

function  results = save_potential_var_parts_to_results(results, mat, num_of_clusters_two_seqs, window, peaks, all_seq, field_names)
% RESULTS = SAVE_POTENTIAL_VAR_PARTS_TO_RESULTS(RESULTS, MAT,
% NUM_OF_CLUSTERS_TWO_SEQS, WINDOW, PEAKS, ALL_SEQ, FIELD_NAMES) Save only results with maximum number of classified clusters.
%   results: structure with variable markers
%   mat: matrix with number of clusters for combination of two variable segments
%   num_of_clusters_two_seqs: vector with number of clusters for combination of two variable segments
%   window: length of words used for entropy calculation
%   peaks: peaks of entropy signals which are above the threshold
%   all_seq: sequences which corresponded to peaks
%   field_names: names of sequences which corresponded to peaks

[row,column] = find(mat == max(num_of_clusters_two_seqs));

for num_of_results = 1:length(row)
    results(end+1).window = window;
    results(end).max_num_of_clusters = max(num_of_clusters_two_seqs);
    results(end).sequence_1 = all_seq.(field_names{row(num_of_results)});
    results(end).positions_1 = peaks(row(num_of_results),:);
    results(end).sequence_2 = all_seq.(field_names{column(num_of_results)});
    results(end).positions_2 = peaks(column(num_of_results),:);
end

function results = get_information_about_segments_from_gb_file(genbank_file, results)
% RESULTS = GET_INFORMATION_ABOUT_SEGMENTS_FROM_GB_FILE(GENBANK_FILE, RESULTS) Get information about variable parts.
%   genbank_file: file name of genbank file with sequence used for reference-based assembly
%   results: structure with variable markers

gb_file = genbankread(genbank_file);

for i = 1:length(gb_file.CDS)
    
    if gb_file.CDS(i).indices(1) < gb_file.CDS(i).indices(2)
        gb_file.CDS(i).start = gb_file.CDS(i).indices(1);
        gb_file.CDS(i).stop = gb_file.CDS(i).indices(2);
    else
        gb_file.CDS(i).stop = gb_file.CDS(i).indices(1);
        gb_file.CDS(i).start = gb_file.CDS(i).indices(2);
    end
end
 
for i = 1:length(results)
    for j = 1:length(gb_file.CDS)
        if gb_file.CDS(j).start <= results(i).positions_1(1) && results(i).positions_1(1) <= gb_file.CDS(j).stop
            results(i).begin_in_gene = gb_file.CDS(j).product;
            results(i).gb_index_start_seq_1 = [gb_file.CDS(j).start gb_file.CDS(j).stop];
        end
        
        if gb_file.CDS(j).start <= results(i).positions_1(2) && results(i).positions_1(2) <= gb_file.CDS(j).stop
            results(i).end_in_gene = gb_file.CDS(j).product;
            results(i).gb_index_end_seq_1 = [gb_file.CDS(j).start gb_file.CDS(j).stop];
        end
        
    end
end

for i = 1:length(results)
    for j = 1:length(gb_file.CDS)
        if gb_file.CDS(j).start <= results(i).positions_2(1) && results(i).positions_2(1) <= gb_file.CDS(j).stop
            results(i).begin_in_gene_2 = gb_file.CDS(j).product;
            results(i).gb_index_start_seq_2 = [gb_file.CDS(j).start gb_file.CDS(j).stop];
        end

        if gb_file.CDS(j).start <= results(i).positions_2(2) && results(i).positions_2(2) <= gb_file.CDS(j).stop
            results(i).end_in_gene_2 = gb_file.CDS(j).product;
            results(i).gb_index_end_seq_2 = [gb_file.CDS(j).start gb_file.CDS(j).stop];
        end

    end
end

function results = get_only_results_with_max_number_of_clusters(results)
% RESULTS = GET_ONLY_RESULTS_WITH_MAX_NUMBER_OF_CLUSTERS(RESULTS) Return only results with maximum number of clusters.
%   results: structure with variable markers

max_num_of_clusters = max([results.max_num_of_clusters]);

i = 1;
while i <= size(results,2)
    
    if results(i).max_num_of_clusters ~= max_num_of_clusters
        results(i) = [];
    else
        i = i + 1;
    end
    
end

