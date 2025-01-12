clear
clc
close all

% filepath: /path/to/your_script.m
% Read the file
filename = 'Thegiantpanda.rtf';
fileID = fopen(filename, 'r');
text = fscanf(fileID, '%c');
fclose(fileID);

% Preprocess the text
text = lower(text); % Convert to lowercase
text = regexprep(text, '[^\w\s]', ''); % Remove punctuation
words = strsplit(text); % Split into words

% Identify unique words and calculate probabilities
wordCounts = containers.Map();
totalWords = length(words);

for i = 1:totalWords
    word = words{i};
    if isKey(wordCounts, word)
        wordCounts(word) = wordCounts(word) + 1;
    else
        wordCounts(word) = 1;
    end
end

wordProbs = containers.Map();
keys = wordCounts.keys;
for i = 1:length(keys)
    word = keys{i};
    wordProbs(word) = wordCounts(word) / totalWords;
end

% Calculate 2-tuple word probabilities
tupleCounts = containers.Map();
for i = 1:totalWords-1
    tuple = [words{i} ' ' words{i+1}];
    if isKey(tupleCounts, tuple)
        tupleCounts(tuple) = tupleCounts(tuple) + 1;
    else
        tupleCounts(tuple) = 1;
    end
end

tupleProbs = containers.Map();
keys = tupleCounts.keys;
for i = 1:length(keys)
    tuple = keys{i};
    tupleProbs(tuple) = tupleCounts(tuple) / (totalWords - 1);
end

% Graphic representation of the 5 most used words
wordCountsArray = cell2mat(values(wordCounts));
[sortedCounts, sortedIndices] = sort(wordCountsArray, 'descend');
top5Words = keys(sortedIndices(1:5));
top5Counts = sortedCounts(1:5);

figure;
bar(categorical(top5Words), top5Counts);
title('Top 5 Most Used Words');
xlabel('Words');
ylabel('Counts');

% Generate a sentence
function sentence = generate_sentence(start_word, wordProbs, tupleProbs, max_length)
    sentence = start_word;
    current_word = start_word;
    for i = 1:max_length-1
        next_word = get_next_word(current_word, tupleProbs);
        if isempty(next_word)
            break;
        end
        sentence = [sentence ' ' next_word];
        current_word = next_word;
    end
end

function next_word = get_next_word(current_word, tupleProbs)
    candidates = {};
    probs = [];
    keys = tupleProbs.keys;
    for i = 1:length(keys)
        tuple = keys{i};
        parts = strsplit(tuple);
        if strcmp(parts{1}, current_word)
            candidates{end+1} = parts{2};
            probs(end+1) = tupleProbs(tuple);
        end
    end
    if isempty(candidates)
        next_word = '';
    else
        % Normalize probabilities
        probs = probs / sum(probs);
        % Generate a random number
        r = rand();
        cumulative_prob = 0;
        for i = 1:length(candidates)
            cumulative_prob = cumulative_prob + probs(i);
            if r <= cumulative_prob
                next_word = candidates{i};
                return;
            end
        end
    end
end

% Example usage
start_word = 'adult';
max_length = 10;
sentence = generate_sentence(start_word, wordProbs, tupleProbs, max_length);
disp(sentence);