%clear all
D=10000; % length of bloom filter
one_n=100; % number of hash functions or alternatively number of ones in an HD vector
dens=one_n/D; % probability of one
N=10000; %total number of elements the initial memory

simul=10; % for Retouched Bloom Filter (RBF)  RBF since we simply simulate it several times


%create HD vectors for each element
HD=zeros(N,D);
for i=1:N
    ind=randperm(D,one_n); %indices of nonzero elements
    HD(i,ind)=1; % this will create sparse binary vectors with exactly one_n ones
end


F_range=50:50:5000; % range of symbols stored in the filters
thr_limit=zeros(1,length(F_range)); %simplifies the life in terms of searching for an optimal threshold 
thr_limit(1,1:30)=20;
thr_limit(1,31:50)=30;
thr_limit(1,51:65)=40;
thr_limit(1,66:80)=50;
thr_limit(1,81:end)=60;

CURRENT_RBF=zeros(length(F_range),3); % store the result for RBF

for ii=1:length(F_range) % for all values in F_range
    disp(ii)
    F=F_range(ii); % set the current number of stored elements
    
    RES1=zeros(simul,1); % collect statistic of TPR for RBF
    RES2=zeros(simul,1); % collect statistic of FPR for RBF
    for s=1:simul
        ind=randperm(N,F);  %random indices to choose elements which will be stored in RBF
        l_ind=zeros(N,1);
        l_ind(ind,1)=1;
        l_ind=find(~l_ind); % indices of elements not stored in RBF
        RBF=double(sum(HD(ind,:),1)>0); % form content of the RBF
        RBF_ind=find(RBF); % find indeces of nonzero elements in RBF
        
        puncture=0.001*(length(RBF_ind)/D); % delete 1% of nonzero elements
        p_ind=randperm(length(RBF_ind),round(D*puncture));  %randomly choose indices to be punctured
        RBF(RBF_ind(p_ind))=0; % delete the chosen indices
        DP=HD*(RBF'); %dot product to all elements
        RES1(s,1)= sum(  DP(ind,1)==one_n )/F; % TPR for RBF
        RES2(s,1)= sum(  DP(l_ind,1)==one_n )/(N-F); %FPR for RBF
    end
    
    CURRENT_RBF(ii,1:2)= mean([RES1,RES2]); % get mean values for TPR and FPR for RBF
    CURRENT_RBF(ii,3)=(CURRENT_RBF(ii,1)+(1-CURRENT_RBF(ii,2)))/2; %  get mean ACC for RBF  
    
    %theoretical estimations of probabilities of elements in the counting filter
    p1=dens;
    %pd_bin = makedist('Binomial','N',F,'p',p1); %create binomial distribution
    bins = 0:F;
    pdf_bins = binopdf(bins,F,p1);
    
    
    %estimate probabiltiy of zero in the thresholded bloom filter
    thr=0:thr_limit(1,ii); %thresholds for setting 1 in the autoscaling bloom filter
    p0_bf=zeros(size(thr));
    for i=1:length(thr)
        p0_bf(1,i)=sum(pdf_bins(1:thr(i)+1));
    end
    
    %Estimate Binomial distributions for given threshold binarization threshold
    bins_2 = 0:one_n;
    
    pdf_bins2=zeros(length(thr), length(bins_2)); %distribution of dot products for true positives
    pdf_bins3=pdf_bins2; %distribution of dot products for false positives
    
    for i=1:length(thr)
        thr_v=thr(i);
        expectation=(F*one_n-  sum(D*([0:i-1]).*pdf_bins(1:i)) )/(F*one_n); % expected density of bloom filter?
        pdf_bins2(i,:) = binopdf(bins_2,one_n,expectation  ); %check distribution for thr 1
        pdf_bins3(i,:) = binopdf(bins_2,one_n,(1-p0_bf(i)));
    end
        
    %calculate TP and FP for a range of parameters: binarizarion threshold (thr) and decision threshold (thr_th)
    TNFP=cell(1,length(thr));
    for i=1:length(thr)
        for thr_th=0:one_n
            
            TNFP{1,i}(thr_th+1,1)= sum(pdf_bins2(i,thr_th+1:end)); % calculate TRP
            TNFP{1,i}(thr_th+1,2)=sum(pdf_bins3(i,thr_th+1:end)); % calculate FRP
            
        end        
        %calculate ACC using TRP and FPR
        TNFP{1,i}(:,3)=(TNFP{1,i}(:,1)+(1-TNFP{1,i}(:,2)))/2;
    end
    
    %brute force search for the best threshold given that TP should be higher than 0.9    
    THR=zeros(length(thr),5); %best threshold for each value
    for i=1:length(thr)
        %disp(i)
        [~,ind]=max(TNFP{1,i}( (TNFP{1,i}(:,1)>=0.9 )  ,3) );
        
        THR(i,1)=ind-1; % get the actual threshold
        THR(i,2:4)=TNFP{1,i}(ind,:); % get the corresponding performance values: TRP, FRP, and ACC
        
    end
    THR(:,5)=thr; % write down the binarization threshold
    
    
    [~,ind]=max(THR(:,4)); % choose the best performing set pf parameters: binarizarion threshld and decision threshold (T) in terms of ACC 
    
    BEST_ABF(ii,:)=THR(ind,:); % write down the performnace of the best performning ABF for the current number of stored elements
    CURRENT_BF(ii,:)=TNFP{1,1}(one_n+1,:); %performance of the nonoptimized BF
    
    
    
    %analytical calculations of performance for BF with optimal parameters (optimized BF) 
    
    one_n_opt=round((D/F)*log(2)  ); % optimal number of hash functions (ones in HD vector)
    if one_n_opt==0 % if it is zero set to one
        one_n_opt=1;
    end
    BEST_BF(ii,1)=one_n_opt; % decision threshold always the chosen one_n_opt
    BEST_BF(ii,2)=1; % TRP is always 1
    BEST_BF(ii,3)=(1- exp(-(one_n_opt *F/D)) )^(one_n_opt); % analytically calculate FPR
    BEST_BF(ii,4)= (1+(1-BEST_BF(ii,3)))/2; %calculate ACC using TRP and FPR
    
    
end


%%

%plot results
figure()

subplot(1,3,1) % plot TPR for all results
hold on
grid on
plot(F_range, BEST_ABF(:,2),'r-.','LineWidth',3) % autoscaling BF
plot(F_range, BEST_BF(:,2),'k','LineWidth',3) % optimized BF
plot(F_range, CURRENT_BF(:,1),'b--','LineWidth',3) % nonoptimized BF
plot(F_range, CURRENT_RBF(:,1),'m:','LineWidth',3) % nonoptimized RBF
% update panel's properties
xlim([0,5000])
ylim([0.88,1])
xlabel('Number of stored elements, \it{n}')
ylabel('TPR')
legend({'autoscaling BF', 'optimized BF', 'nonoptimized BF', 'optimized RBF'})
subplot(1,3,2) % plot FPR for all results
hold on
grid on
plot(F_range, BEST_ABF(:,3),'r-.','LineWidth',3) % autoscaling BF
plot(F_range, BEST_BF(:,3),'k','LineWidth',3) % optimized BF
plot(F_range, CURRENT_BF(:,2),'b--','LineWidth',3) % nonoptimized BF
plot(F_range, CURRENT_RBF(:,2),'m:','LineWidth',3) % nonoptimized RBF
% update panel's properties
xlim([0,5000])
xlabel('Number of stored elements, \it{n}')
ylabel('FPR')
subplot(1,3,3)  %plot ACC for all results
hold on
grid on
plot(F_range, BEST_ABF(:,4),'r-.','LineWidth',3) % autoscaling BF
plot(F_range, BEST_BF(:,4),'k','LineWidth',3) % optimized BF
plot(F_range, CURRENT_BF(:,3),'b--','LineWidth',3) % nonoptimized BF
plot(F_range, CURRENT_RBF(:,3),'m:','LineWidth',3) % nonoptimized RBF
% update panel's properties
xlim([0,5000])
xlabel('Number of stored elements, \it{n}')
ylabel('ACC')
ylim([0.49,1])

