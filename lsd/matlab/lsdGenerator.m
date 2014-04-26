%
% Input by users (Warning, for best results only use primers/tags of same
% taxa, to avoid the algorithm being skewed by fake inter-tags distances)
%
nbSamples = 100;
% Fasta file for the primers
primerFastaFile = 'primersForams.fasta';
% Wildcards to match with forward primers
forwardPrimerWildcard = {'V9F-', 'F1-', 'V4F-'};
% Wildcards to match with reverse primers
reversePrimerWildcard = {'sBnew-', '15-', 'V4R-'};
%
% Analysis part
%
% Import the corresponding tagged primers file
primersInfo = fastaread(primerFastaFile);
nbForwardTags = 0;
forwardSequences = cell(length(primersInfo), 1);
forwardHeaders = cell(length(primersInfo), 1);
nbReverseTags = 0;
reverseSequences = cell(length(primersInfo), 1);
reverseHeaders = cell(length(primersInfo), 1);
% Perform the analysis
for p = 1:length(primersInfo)
    curHeader = primersInfo(p).Header;
    curSequence = primersInfo(p).Sequence;
    forwardMatch = sum(cell2mat(regexp(curHeader, forwardPrimerWildcard, 'once')));
    reverseMatch = sum(cell2mat(regexp(curHeader, reversePrimerWildcard, 'once')));
    if (forwardMatch && reverseMatch)
        disp(['Error : Found primer "' curHeader '" is a match to both forward and reverse wildcards']);
        return;
    end
    if (forwardMatch)
        nbForwardTags = nbForwardTags + 1;
        forwardSequences{nbForwardTags} = curSequence;
        forwardHeaders{nbForwardTags} = curHeader;
    end
    if (reverseMatch)
        nbReverseTags = nbReverseTags + 1;
        reverseSequences{nbReverseTags} = curSequence;
        reverseHeaders{nbReverseTags} = curHeader;
    end
end
% Discard the extra empty sapce
forwardSequences = forwardSequences(1:nbForwardTags);
reverseSequences = reverseSequences(1:nbReverseTags);
% Analyze the per-sequence edit distances of forward tags
forwardDistances = zeros(nbForwardTags);
for i = 1:nbForwardTags
    for j = (i + 1):nbForwardTags
        score = levenshtein(forwardSequences{i}, forwardSequences{j});
        forwardDistances(i, j) = score;
        forwardDistances(j, i) = score;
        if (score < 3)
            disp('Warning : Two forward primers found at less than 3 edit differences :');
            disp(forwardSequences{i});
            disp(forwardSequences{j});
        end
    end
end
totalForwardDistances = sum(forwardDistances);
% Analyze the per-sequence edit distances of reverse tags
reverseDistances = zeros(nbReverseTags);
for i = 1:nbReverseTags
    for j = (i + 1):nbReverseTags
        score = levenshtein(reverseSequences{i}, reverseSequences{j});
        reverseDistances(i, j) = score;
        reverseDistances(j, i) = score;
        if (score < 3)
            disp('Warning : Two reverse primers found at less than 3 edit differences :');
            disp(reverseSequences{i});
            disp(reverseSequences{j});
        end
    end
end
totalReverseDistances = sum(forwardDistances);
disp(nbForwardTags);
disp(nbReverseTags);
% First find the matrix order
order = lsdCombinatorial(min(nbForwardTags, nbReverseTags), 1);
% If there are more forward than reverse tags, we select the most distinct
if (order < nbForwardTags)
    [totalForwardDistances, distID] = sort(totalForwardDistances, 'descend');
    finalIDs = distID(1:order);
    totalForwardDistances = totalForwardDistances(finalIDs);
    forwardSequences = forwardSequences(finalIDs);
    forwardHeaders = forwardHeaders(finalIDs);
end
% If there are more forward than reverse tags, we select the most distinct
if (order < nbReverseTags)
	[totalReverseDistances, distID] = sort(totalReverseDistances, 'descend');
    finalIDs = distID(1:order);
    reverseSequences = reverseSequences(finalIDs);
    reverseHeaders = reverseHeaders(finalIDs);
end
% Number of potential combinations
nbCombinations = order * order;
% Number of samples per tags
samplePerTags = ceil(nbSamples / order);
if (samplePerTags <= 1)
    disp('Your amount of tags allow for a fully non-combinatorial design');
    return;
end
% Warn user about saturation effect if more than half of samples are used
if (nbSamples > (nbCombinations / 2)), disp('Warning : Potential saturation effect'); end
if (nbSamples > ((3 * nbCombinations) / 4)), disp('Warning : Almost certain saturation'); end
if (nbSamples == nbCombinations), disp('Warning : Complete saturation, mistagging impossible to detect'); end
% First retrieve the complete set of Latin Square Designs from specific order
disp('Generating candidate matrix.');
disp(order);
candidateMatrix = lsdCombinatorial(order, 0);
disp(['  = Computed ' num2str(length(candidateMatrix)) ' potential matrix arrangements']);
if (length(candidateMatrix) > 10000)
    candidateMatrix = datasample(candidateMatrix, 10000, 'Replace', false);
    disp('  = Matrix pool too large, limiting to random set of 10000');
end
% Now we need to select symbols in each potential solutions
disp('Computing symbol permutations.');
symbolPermutations = combnk(1:order, samplePerTags);
% Number of potential solutions
nbSolutions = length(candidateMatrix) * size(symbolPermutations, 1);
disp(['  = Found ' num2str(nbSolutions) ' potential solutions']);
if (nbSolutions > 1000000)
    disp('  = Solutions pool too large, limiting to random 1000000');
    nbMaxPermutations = ceil(100000 / length(candidateMatrix));
    symbolPermutations = datasample(symbolPermutations, nbMaxPermutations, 1, 'Replace', false);
    nbSolutions = length(candidateMatrix) * size(symbolPermutations, 1);
end
% We are going to generate every potential solution matrix
solutionMatrix = cell(nbSolutions, 1);
% We want to keep track of their fitness (how well they fit the constraints)
solutionFitness = zeros(nbSolutions, 1);
% Current analyzed solution
curSolution = 1;
disp('Evaluating fitness of solutions. (Dot every 10% chunks)');
% Now parse through every potential matrix
for m = 1:length(candidateMatrix)
	curMatrix = candidateMatrix{m};
    if (mod(10 * m / length(candidateMatrix), 1) == 0)
        fprintf('.');
    end
	% Analyze every potential permutation of symbols
	for p = 1:size(symbolPermutations, 1)
		curPerm = symbolPermutations(p, :);
		finalMatrix = zeros(order, order);
		for e = 1:length(curPerm)
			finalMatrix(curMatrix == curPerm(e)) = 1;
		end
		solutionMatrix{curSolution} = finalMatrix;
		solutionFitness(curSolution) = lsdFitness(finalMatrix);
		curSolution = curSolution + 1;
	end
end
fprintf('\n');
% Now sort solutions by fitness
[fitness solIDs] = sort(solutionFitness, 'descend');
% We generated a set of solutions based on a purely equal usage frequency
nbSamplesGenerated = samplePerTags * order;
% Number of samples to remove from the best solutions
nbExcessSamples = nbSamplesGenerated - nbSamples;
disp(['Removing the excess of ' num2str(nbExcessSamples) ' samples.']);
% Now we only process the 10 best solutions
for i = 1:10
    matrix = solutionMatrix{solIDs(i)};
    % We will start by removing the "denser samples"
    [sampleX, sampleY] = find(matrix);
    sampleCoords = [sampleX sampleY];
    distMatrix = pdist(sampleCoords, 'euclidean');
    fullDistMatrix = sum(squareform(distMatrix));
    for s = nbExcessSamples:-1:1
        [~, xPos] = min(fullDistMatrix);
        matrix(sampleCoords(xPos, 1), sampleCoords(xPos, 2)) = 0;
        fullDistMatrix(xPos) = Inf;
    end
    % Display the best matrix
    disp(matrix);
    disp(sum(sum(matrix)));
end