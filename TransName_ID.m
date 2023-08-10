clear
model1 = readCbModel('iEC1356_Bl21DE3','fileType','Matlab');
model=model1;
% model = readCbModel('iEC1356_Bl21DE3','fileType','SBML');
%% Read the FASTA file to get the ID and sequence information
[id1,seq1] = fastaread('K5T46.fna');  
[id2,seq2] = fastaread('ECD_RS.fna');
[id3,seq3] = fastaread('ECD_000.fna'); 
%% Segmentation of gene annotation information   
for i=1:length(id1)
    s=id1{i};
first_idx = strfind(s, '[locus_tag=');
last_idx = strfind(s, '] [p');
locus_tag1{i,1} = (s(first_idx:last_idx(1)));

first_idx = strfind(s, '[protein='); 
last_idx = strfind(s, '] [p');
protein1{i,1} = (s(first_idx:last_idx(2)));

first_idx = strfind(s, '[protein_id=');  
last_idx = strfind(s, '] [location'); 
protein_id1{i,1} = (s(first_idx:last_idx));
end


for i=1:length(id2)
    s=id2{i};
first_idx = strfind(s, '[locus_tag=');
last_idx = strfind(s, '] [p');
locus_tag2{i,1} = (s(first_idx:last_idx(1)));

first_idx = strfind(s, '[protein='); 
last_idx = strfind(s, '] [p');
protein2{i,1} = (s(first_idx:last_idx(2)));

first_idx = strfind(s, '[protein_id=');  
last_idx = strfind(s, '] [location'); 
protein_id2{i,1} = (s(first_idx:last_idx));
end


for i=1:length(locus_tag1) 
str = locus_tag1{i};
start_idx = regexp(str, '\[locus_tag=')   ;
end_idx = regexp(str, '\]')         ;     
match= str(start_idx+11:end_idx-1);
tag1{i,1}=match;
end

for i=1:length(locus_tag2) 
str = locus_tag2{i};
start_idx = regexp(str, '\[locus_tag=')   ;
end_idx = regexp(str, '\]')         ;     
match= str(start_idx+11:end_idx-1);
tag2{i,1}=match;
end

%% 匹配序列
seq1=seq1';
F1=cell(length(seq1), 4);
for i=1:length(seq1)   
    F1{i,2}= tag2((ismember(seq2,seq1{i,1})),1);
    F1{i,1} = tag1(i,1);    
    F1{i,3} = protein1{i};  
    F1{i,4} = protein_id1{i};  
end
T1 = table(F1(:,1),F1(:,2),F1(:,3),F1(:,4));
T1.Properties.VariableNames = {'locus_tag1', 'locus_tag2','protein', 'protein_id'};


seq2=seq2';

F2=cell(length(seq2), 4);
for i=1:length(seq2)    
    F2{i,1}  = tag1((ismember(seq1,seq2{i,1})),1);    
    F2{i,2} = tag2(i,1); 
    F2{i,3} = protein2{i};  
    F2{i,4} = protein_id2{i};   
end

T2 = table(F2(:,1),F2(:,2),F2(:,3),F2(:,4));
T2.Properties.VariableNames = {'locus_tag1', 'locus_tag2','protein', 'protein_id'};

%% map到模型中的基因
Search2={};
for i=1:length(model.genes)
Search2{i,1}=  model.genes{i,1};
Temp2= F2((ismember(tag2,model.genes{i,1})),1);
if isempty(Temp2)==1
    Search2{i,2}='nan';
else
    Search2{i,2}=Temp2{1,1};
end
Temp3=F2((ismember(tag2,model.genes{i,1})),4);
if isempty(Temp3)==1
    Search2{i,3}='nan';
else
    Search2{i,3}=Temp3{1,1};
end
end

