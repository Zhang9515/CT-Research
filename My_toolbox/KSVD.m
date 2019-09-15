function [Dictionary,output] = KSVD(...  

   Data,... % an nXN matrix that contins N signals (Y), each of dimension n.  

   param)  

% =========================================================================  

%                          K-SVD algorithm  

% =========================================================================  

% The K-SVD algorithm finds a dictionary for linear representation of  

% signals. Given a set of signals, it searches for the best dictionary that  

% can sparsely represent each signal. Detailed discussion on the algorithm  

% and possible applications can be found in "The K-SVD: An Algorithm for  

% Designing of Overcomplete Dictionaries for Sparse Representation", written  

% by M. Aharon, M. Elad, and A.M. Bruckstein and appeared in the IEEE Trans.  

% On Signal Processing, Vol. 54, no. 11, pp. 4311-4322, November 2006.  

% =========================================================================  

% INPUT ARGUMENTS:  

% Data                         an nXN matrix that contins N signals (Y), each of dimension n.  

% param                        structure that includes all required  

%                                 parameters for the K-SVD execution.  

%                                 Required fields are:  

%    K, ...                    the number of dictionary elements to train  

%    numIteration,...          number of iterations to perform.  

%    errorFlag...              if =0, a fix number of coefficients is  

%                                 used for representation of each signal. If so, param.L must be  

%                                 specified as the number of representing atom. if =1, arbitrary number  

%                                 of atoms represent each signal, until a specific representation error  

%                                 is reached. If so, param.errorGoal must be specified as the allowed  

%                                 error.  

%    preserveDCAtom...         if =1 then the first atom in the dictionary  

%                                 is set to be constant, and does not ever change. This  

%                                 might be useful for working with natural  

%                                 images (in this case, only param.K-1  

%                                 atoms are trained).  

%    (optional, see errorFlag) L,...                 % maximum coefficients to use in OMP coefficient calculations.  

%    (optional, see errorFlag) errorGoal, ...        % allowed representation error in representing each signal.  

%    InitializationMethod,...  mehtod to initialize the dictionary, can  

%                                 be one of the following arguments:  

%                                 * 'DataElements' (initialization by the signals themselves), or:  

%                                 * 'GivenMatrix' (initialization by a given matrix param.initialDictionary).  

%    (optional, see InitializationMethod) initialDictionary,...      % if the initialization method  

%                                 is 'GivenMatrix', this is the matrix that will be used.  

%    (optional) TrueDictionary, ...        % if specified, in each  

%                                 iteration the difference between this dictionary and the trained one  

%                                 is measured and displayed.  

%    displayProgress, ...      if =1 progress information is displyed. If param.errorFlag==0,  

%                                 the average repersentation error (RMSE) is displayed, while if  

%                                 param.errorFlag==1, the average number of required coefficients for  

%                                 representation of each signal is displayed.  

% =========================================================================  

% OUTPUT ARGUMENTS:  

%  Dictionary                  The extracted dictionary of size nX(param.K).  

%  output                      Struct that contains information about the current run. It may include the following fields:  

%    CoefMatrix                  The final coefficients matrix (it should hold that Data equals approximately Dictionary*output.CoefMatrix.  

%    ratio                       If the true dictionary was defined (in  

%                                synthetic experiments), this parameter holds a vector of length  

%                                param.numIteration that includes the detection ratios in each  

%                                iteration).  

%    totalerr                    The total representation error after each  

%                                iteration (defined only if  

%                                param.displayProgress=1 and  

%                                param.errorFlag = 0)  

%    numCoef                     A vector of length param.numIteration that  

%                                include the average number of coefficients required for representation  

%                                of each signal (in each iteration) (defined only if  

%                                param.displayProgress=1 and  

%                                param.errorFlag = 1)  

% =========================================================================  

 

if (~isfield(param,'displayProgress'))  

   param.displayProgress = 0;  

end  

totalerr(1) = 99999;  

if (isfield(param,'errorFlag')==0)  

   param.errorFlag = 0;  

end  

 

if (isfield(param,'TrueDictionary'))  

   displayErrorWithTrueDictionary = 1;  

   ErrorBetweenDictionaries = zeros(param.numIteration+1,1);  

   ratio = zeros(param.numIteration+1,1);  

else  

   displayErrorWithTrueDictionary = 0;  

   ratio = 0;  

end  

if (param.preserveDCAtom>0)  

   FixedDictionaryElement(1:size(Data,1),1) = 1/sqrt(size(Data,1));  

else  

   FixedDictionaryElement = [];  

end  

% coefficient calculation method is OMP with fixed number of coefficients  

 

if (size(Data,2) < param.K)  

   disp('Size of data is smaller than the dictionary size. Trivial solution...');  

   Dictionary = Data(:,1:size(Data,2));  

   return;  

   %若参数K大于信号的个数 则将数据集作为字典集  

elseif (strcmp(param.InitializationMethod,'DataElements'))  

   Dictionary(:,1:param.K-param.preserveDCAtom) = Data(:,1:param.K-param.preserveDCAtom);%将数据集的1到param.K-param.preserveDCAtom数据作为字典集  

elseif (strcmp(param.InitializationMethod,'GivenMatrix'))  

   Dictionary(:,1:param.K-param.preserveDCAtom) = param.initialDictionary(:,1:param.K-param.preserveDCAtom);%将initialDictionary的1到param.K-param.preserveDCAtom列作为字典集  

end  

% reduce the components in Dictionary that are spanned by the fixed  

% elements  

if (param.preserveDCAtom)  

   tmpMat = FixedDictionaryElement Dictionary;%求minimal norm（FixedDictionaryElement×Dictionary-dictionary）  

   Dictionary = Dictionary - FixedDictionaryElement*tmpMat;  

end  

%normalize the dictionary.  

Dictionary = Dictionary*diag(1./sqrt(sum(Dictionary.*Dictionary)));%归一化  

Dictionary = Dictionary.*repmat(sign(Dictionary(1,:)),size(Dictionary,1),1); % multiply in the sign of the first element.%字典集中的每个元素的化为正数  

totalErr = zeros(1,param.numIteration);  

 

% the K-SVD algorithm starts here.  

 

for iterNum = 1:param.numIteration  

   % find the coefficients  

   if (param.errorFlag==0)%固定表达系数的个数  

       %CoefMatrix = mexOMPIterative2(Data, [FixedDictionaryElement,Dictionary],param.L);  

       CoefMatrix = OMP([FixedDictionaryElement,Dictionary],Data, param.L);  

   else  

       %CoefMatrix = mexOMPerrIterative(Data, [FixedDictionaryElement,Dictionary],param.errorGoal);  

       CoefMatrix = OMPerr([FixedDictionaryElement,Dictionary],Data, param.errorGoal);%设定允许的误差  

       param.L = 1;  

   end  

     

   replacedVectorCounter = 0;  

   rPerm = randperm(size(Dictionary,2));%生成一个1到size(Dictionary,2)的随机的向量  

   for j = rPerm  

       [betterDictionaryElement,CoefMatrix,addedNewVector] = I_findBetterDictionaryElement(Data,...  

           [FixedDictionaryElement,Dictionary],j+size(FixedDictionaryElement,2),...  

           CoefMatrix ,param.L);  

       Dictionary(:,j) = betterDictionaryElement;%更新字典集  

       if (param.preserveDCAtom)  

           tmpCoef = FixedDictionaryElementbetterDictionaryElement;  

           Dictionary(:,j) = betterDictionaryElement - FixedDictionaryElement*tmpCoef;  

           Dictionary(:,j) = Dictionary(:,j)./sqrt(Dictionary(:,j)'*Dictionary(:,j));%归一化  

       end  

       replacedVectorCounter = replacedVectorCounter+addedNewVector;  

   end  

 

   if (iterNum>1 & param.displayProgress)  

       if (param.errorFlag==0)  

           output.totalerr(iterNum-1) = sqrt(sum(sum((Data-[FixedDictionaryElement,Dictionary]*CoefMatrix).^2))/prod(size(Data)));  

           disp(['Iteration   ',num2str(iterNum),'   Total error is: ',num2str(output.totalerr(iterNum-1))]);  

       else  

           output.numCoef(iterNum-1) = length(find(CoefMatrix))/size(Data,2);  

           disp(['Iteration   ',num2str(iterNum),'   Average number of coefficients: ',num2str(output.numCoef(iterNum-1))]);  

       end  

   end  

   if (displayErrorWithTrueDictionary )  

       [ratio(iterNum+1),ErrorBetweenDictionaries(iterNum+1)] = I_findDistanseBetweenDictionaries(param.TrueDictionary,Dictionary);  

       disp(strcat(['Iteration  ', num2str(iterNum),' ratio of restored elements: ',num2str(ratio(iterNum+1))]));  

       output.ratio = ratio;  

   end  

   Dictionary = I_clearDictionary(Dictionary,CoefMatrix(size(FixedDictionaryElement,2)+1:end,:),Data);  

     

   if (isfield(param,'waitBarHandle'))  

       waitbar(iterNum/param.counterForWaitBar);  

   end  

end  

 

output.CoefMatrix = CoefMatrix;  

Dictionary = [FixedDictionaryElement,Dictionary];  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%  findBetterDictionaryElement  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

 

function [betterDictionaryElement,CoefMatrix,NewVectorAdded] = I_findBetterDictionaryElement(Data,Dictionary,j,CoefMatrix,numCoefUsed)%CoefMatrix为字典最终的系数  

if (length(who('numCoefUsed'))==0)  

   numCoefUsed = 1;  

end  

relevantDataIndices = find(CoefMatrix(j,:)); % 非零元在第j行的系数矩阵中的位置the data indices that uses the j'th dictionary element.  

if (length(relevantDataIndices)<1) %(length(relevantDataIndices)==0)如果系数矩阵的第j列全为零  

   ErrorMat = Data-Dictionary*CoefMatrix;%在已有的字典集下和系数下对data项的估计误差  

   ErrorNormVec = sum(ErrorMat.^2);%对误 差每项平方  

   [d,i] = max(ErrorNormVec);%d为所有列中最大项，i为其第几列  

   betterDictionaryElement = Data(:,i);%ErrorMat(:,i); %数据项的i列赋给betterDictionaryElement  

   betterDictionaryElement = betterDictionaryElement./sqrt(betterDictionaryElement'*betterDictionaryElement);%归一化betterDictionaryElement  

   betterDictionaryElement = betterDictionaryElement.*sign(betterDictionaryElement(1));%将betterDictionaryElement中负的元素化为正的  

   CoefMatrix(j,:) = 0;%将系数矩阵的第j行赋值为0  

   NewVectorAdded = 1;  

   return;  

end  

 

NewVectorAdded = 0;  

tmpCoefMatrix = CoefMatrix(:,relevantDataIndices); %将系数矩阵的第j行的非零项所在的列赋给tmpCoefMatrix  

tmpCoefMatrix(j,:) = 0;% the coeffitients of the element we now improve are not relevant.将tmpCoefMatrix第j行赋0  

errors =(Data(:,relevantDataIndices) - Dictionary*tmpCoefMatrix); %在除去字典中第j个的元素后数据集与预测数据之间的误差% vector of errors that we want to minimize with the new element  

% % the better dictionary element and the values of beta are found using svd.  

% % This is because we would like to minimize || errors - beta*element ||_F^2.  

% % that is, to approximate the matrix 'errors' with a one-rank matrix. This  

% % is done using the largest singular value.  

[betterDictionaryElement,singularValue,betaVector] = svds(errors,1);%betterDictionaryElement为右特征向量 singularValue为最大特征值 betaVector左特征向量  

CoefMatrix(j,relevantDataIndices) = singularValue*betaVector';% *signOfFirstElem 系数矩阵的第j行的非零元的位置换为singularValue*betaVector的值  

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%  findDistanseBetweenDictionaries  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

function [ratio,totalDistances] = I_findDistanseBetweenDictionaries(original,new)  

% first, all the column in oiginal starts with positive values.  

catchCounter = 0;  

totalDistances = 0;  

for i = 1:size(new,2)  

   new(:,i) = sign(new(1,i))*new(:,i);  

end  

for i = 1:size(original,2)  

   d = sign(original(1,i))*original(:,i);  

   distances =sum ( (new-repmat(d,1,size(new,2))).^2);  

   [minValue,index] = min(distances);  

   errorOfElement = 1-abs(new(:,index)'*d);  

   totalDistances = totalDistances+errorOfElement;  

   catchCounter = catchCounter+(errorOfElement<0.01);  

end  

ratio = 100*catchCounter/size(original,2);  

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%  I_clearDictionary  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

function Dictionary = I_clearDictionary(Dictionary,CoefMatrix,Data)  

T2 = 0.99;  

T1 = 3;  

K=size(Dictionary,2);  

Er=sum((Data-Dictionary*CoefMatrix).^2,1); % remove identical atoms  

G=Dictionary'*Dictionary; G = G-diag(diag(G));  

for jj=1:1:K,  

   if max(G(jj,:))>T2 | length(find(abs(CoefMatrix(jj,:))>1e-7))<=T1 ,  

       [val,pos]=max(Er);  

       Er(pos(1))=0;  

       Dictionary(:,jj)=Data(:,pos(1))/norm(Data(:,pos(1)));  

       G=Dictionary'*Dictionary; G = G-diag(diag(G));  

   end;  

end;  