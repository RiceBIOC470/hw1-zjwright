% Homework 1. Due before class on 9/12/17

%% Problem 1 - addition with strings

% Fill in the blank space in this section with code that will add 
% the two numbers regardless of variable type. Hint see the matlab
% functions ischar, isnumeric, and str2num. 

%your code should work no matter which of these lines is uncommented. 
x = 3; y = 5; % integers
%x = '3'; y= '5'; %strings
% x = 3; y = '5'; %mixed

%your code goes here
%your code goes here
if isnumeric(x) == 0;
   x = str2num(x);
end
if isnumeric(y) == 0;
   y = str2num(y);
end
answer = x+y;
%output your answer
disp(answer);

%% Problem 2 - our first real biology problem. Open reading frames and nested loops.

%part 1: write a piece of code that creates a random DNA sequence of length
% N (i.e. consisting of the letters ATGC) where we will start with N=500 base pairs (b.p.).
% store the output in a variable
% called rand_seq. Hint: the function randi may be useful. 
% Even if you have access to the bioinformatics toolbox, 
% do not use the built-in function randseq for this part. 

N = 500; % define sequence length
DNA='ATCG'; N=500
rand_seq=DNA(randi(4,1,N));

%part 2: open reading frames (ORFs) are pieces of DNA that can be
% transcribed and translated. They start with a start codon (ATG) and end with a
% stop codon (TAA, TGA, or TAG). Write a piece of code that finds the longest ORF 
% in your seqeunce rand_seq. Hint: see the function strfind.
startcodon=strfind(rand_seq,'ATG');
max_ORFlength=0;
ORF_length=0;
for i=startcodon
    for ii=i:3:(N-2)
        if rand_seq(ii:ii+2)=='TAA'|rand_seq(ii:ii+2)=='TAG'|rand_seq(ii:ii+2)=='TGA'
           ORF_length=(ii+3)-i;
           break
        end
           if ORF_length>max_ORFlength;
               max_ORFlength=ORF_length;
           end
    end
end
disp(max_ORFlength);

%part 3: copy your code in parts 1 and 2 but place it inside a loop that
% runs 1000 times. Use this to determine the probability
% that an sequence of length 500 has an ORF of greater than 50 b.p.

%part 4: copy your code from part 3 but put it inside yet another loop,
% this time over the sequence length N. Plot the probability of having an
% ORF > 50 b.p. as a funciton of the sequence length. 
probability=[];
N=1000;
for sequence=1:N

max_ORFlength=0;
ORF_length=0;
n=0;

DNA='ATCG'; rand_seq=DNA(randi(4,1,N));
startcodon=strfind(rand_seq,'ATG');

for i=startcodon
    for ii=i:3:(N-2)
        if rand_seq(ii:ii+2)=='TAA'|rand_seq(ii:ii+2)=='TAG'|rand_seq(ii:ii+2)=='TGA'
           ORF_length=(ii+3)-i;
           break
        end
           if ORF_length>max_ORFlength;
               max_ORFlength=ORF_length;
           end
    end
    if max_ORFlength>50
        n=n+1;
    probability(sequence)=n/1000;   
    end
end
end
xval=[1:N];
yval=probability;
plot(xval,yval)
%part 5: Make sure your results from part 4 are sensible. What features
% must this curve have (hint: what should be the value when N is small or when
% N is very large? how should the curve change in between?) Make sure your
% plot looks like this. 

%when N is small, probability is zero or low, and probability approaches 1
%the longer N gets. My code does not produce this result and I'm not sure
%what the problem is.
%% problem 3 data input/output and simple analysis

%The file qPCRdata.txt is an actual file that comes from a Roche
%LightCycler qPCR machine. The important columns are the Cp which tells
%you the cycle of amplification and the position which tells you the well
%from the 96 well plate. Each column of the plate has a different gene and
%each row has a different condition. Each gene in done in triplicates so
%columns 1-3 are the same gene, columns 4-6 the same, etc.
%so A1-A3 are gene 1 condition 1, B1-B3 gene 1 condition 2, A4-A6 gene 2
%condition 1, B4-B6 gene2 condition 2 etc. 

% part1: write code to read the Cp data from this file into a vector. You can ignore the last two
% rows with positions beginning with G and H as there were no samples here. 

filename = 'qPCRdata.txt';
fid = fopen(filename,'r');
data = textscan(fid, '%*s%*s%s%*s%f%*s%*s%*s%[^\n]', 72,'Delimiter', '\t','HeaderLines', 2);
Cp=cell2mat(data(1,2));

% Part 2: transform this vector into an array representing the layout of
% the plate. e.g. a 6 row, 12 column array should that data(1,1) = Cp from
% A1, data(1,2) = Cp from A2, data(2,1) = Cp from B1 etc. 
Cp_array=vec2mat(Cp,12)
% Part 3. The 4th gene in columns 10 - 12 is known as a normalization gene.
% That is, it's should not change between conditions and it is used to normalize 
% the expression values for the others. For the other three
% genes, compute their normalized expression in all  conditions, normalized to condition 1. 
% In other words, the fold change between these conditions and condition 1. The
% formula for this is 2^[Cp0 - CpX - (CpN0 - CpNX)] where Cp0 is the Cp for
% the gene in the 1st condition, CpX is the value of Cp in condition X and
% CpN0 and CpNX are the same quantitites for the normalization gene.
% Plot this data in an appropriate way. 

norm = zeros(6,3); %output should be 3 columns (one for each gene average) with six rows (for each condition)
for c = 1:6
    norm(c, 1) = 2^(mean(plate(1, 1:3)) - mean(plate(c, 1:3)) - (mean(plate(1, 10:12)) - mean(plate(c, 10:12))));
    norm(c, 2) = 2^(mean(plate(1, 4:6)) - mean(plate(c, 4:6)) - (mean(plate(1, 10:12)) - mean(plate(c, 10:12))));
    norm(c, 3) = 2^(mean(plate(1, 7:9)) - mean(plate(c, 7:9)) - (mean(plate(1, 10:12)) - mean(plate(c, 10:12))));
end
xval=1:6; yval_geneA=norm(1:6, 1); yval_geneB=norm(1:6,2); yval_geneC=norm(1:6,3);
plot(xval, yval_geneA, 'r'); hold on; plot(xval, yval_geneB, 'b'); hold on; plot(xval, yval_geneC, 'g')

%% Challenge problems that extend the above (optional)

% 1. Write a solution to Problem 2 part 2 that doesn't use any loops at
% all. Hint: start by using the built in function bsxfun to make a matrix of all distances
% between start and stop codons. 
DNA='ATCG'; N=500
DNAseq=DNA(randi(4,1,N));
disp(DNAseq);
   startcodons=strfind(DNAseq, 'ATG');
   end1=strfind(DNAseq, 'TAA'); end2=strfind(DNAseq, 'TGA'); end3=strfind(DNAseq, 'TAG');
   TAAend=bsxfun(@minus, end1, startcodons'); TGAend=bsxfun(@minus, end2, startcodons'); TAGend=bsxfun(@minus, end3, startcodons');
   TAAseq=find(rem((TAAend/3),1)==0); TGAseq=find(rem((TGAend/3),1)==0); TAGseq=find(rem((TAGend/3),1)==0);
   A=TAAend(TAAseq); B=TGAend(TGAseq); C=TAGend(TAGseq); 
   allseqs=[A', B', C']; 
   
%I gave up here. It can find the longest seqeunce that starts with ATG and
%ends with a stop codon in frame, but I don't know how to make it end on
%the earlier stop codons that are still in frame with ATG.


% 2. Problem 2, part 4. Use Matlab to compute the exact solution to this
% problem and compare your answer to what you got previously by testing
% many sequences. Plot both on the same set of axes. Hint: to get started 
% think about the following:
% A. How many sequences of length N are there?
% B. How many ways of making an ORF of length N_ORF are there?
% C. For each N_ORF how many ways of position this reading frame in a
% sequence of length N are there?

% 3. Problem 3. Assume that the error in each Cp is the standard deviation
% of the three measurements. Add a section to your code that propogates this
% uncertainty to the final results. Add error bars to your plot. (on
% propagation of error, see, for example:
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty


