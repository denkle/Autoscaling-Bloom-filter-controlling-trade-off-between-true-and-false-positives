clear all
D=10000; % length of bloom filter 
one_n=100; % number of hash functions or alternatively number of ones in an HD vector
dens=one_n/D; % probability of one
F=500; %total number of elements stored the filter
%N=10000; %total number of elements the initial memory
%disp((D/F)*log(2)) % optimal number of hash functions in bloom filter for the given parameters


%theoretical estimations of probabilities of elements in the counting filter
p1=dens;
%pd_bin = makedist('Binomial','N',F,'p',p1); %create binomial distribution
bins = 0:F;
pdf_bins = binopdf(bins,F,p1); 


%estimate probabiltiy of zero in the thresholded bloom filter
thr=0:5; %thresholds for setting 1 in the autoscaling bloom filter
p0_bf=zeros(size(thr));
for i=1:length(thr)
    p0_bf(1,i)=sum(pdf_bins(1:thr(i)+1));
end

%Estimate Binomial distributions for given binarization threshold
bins_2 = 0:one_n;
for i=1:length(thr)
    thr_v=thr(i);
    expectation=(F*one_n  - sum(D*([0:i-1]).*pdf_bins(1:i))      )/(F*one_n); % expected density of bloom filter?
    pdf_bins2(i,:) = binopdf(bins_2,one_n,expectation  ); %distribution of dot products for true positives
    pdf_bins3(i,:) = binopdf(bins_2,one_n,(1-p0_bf(i))); %distribution of dot products for false positives
end

%calculate TP and FP for a range of parameters: binarizarion threshold
%(thr) and decision threshold (thr_th)
TNFP=cell(1,length(thr));
for i=1:length(thr)
    for thr_th=0:one_n
        
        TNFP{1,i}(thr_th+1,1)= sum(pdf_bins2(i,thr_th+1:end)); % calculate TRP
        TNFP{1,i}(thr_th+1,2)=sum(pdf_bins3(i,thr_th+1:end)); %calculate FRP
        
    end
    
    %calculate ACC using TRP and FPR
    TNFP{1,i}(:,3)=(TNFP{1,i}(:,1)+(1-TNFP{1,i}(:,2)))/2;
end

%brute force search for the best threshold given that TP should be higher than 0.97
THR=zeros(length(thr),4); %best threshold for each value
for i=1:length(thr)
    [~,ind]=max(TNFP{1,i}( (TNFP{1,i}(:,1)>=0.97 )  ,3) );
    
    THR(i,1)=ind-1;% get the actual threshold
    THR(i,2:4)=TNFP{1,i}(ind,:); % get the corresponding performance values: TRP, FRP, and ACC
    
end

%plot results
figure()
for i=1:length(thr)
    subplot(2,3,i)
    box on
    hold on
    stem(bins_2,pdf_bins2(i,:),'b','LineWidth',2) % plot distribution of dot products for true positives 
    stem(bins_2,pdf_bins3(i,:),'r','LineWidth',2) % plot distribution of dot products for false positives
    ampl=0.5* max( [pdf_bins2(i,:),pdf_bins3(i,:)] );% plot amplitude bar for threshold T
    plot([THR(i,1),THR(i,1)], [0,ampl],'k','LineWidth',3)
    xlabel('Dot product')
    ylabel('Probability')
    
    %used to set axis limits and legends
    if i==1
        xlim([95,100])
        ylim([0.0,1.0])
    elseif i==2
        xlim([90,100])
        ylim([0.0,0.6])
    elseif i==3
        xlim([80,100])
        ylim([0.0,0.2])
    elseif i==4
        xlim([40,100])
        ylim([0.0,0.15])
    elseif i==5
        xlim([40,100])
        ylim([0.0,0.1])
    elseif i==6
        xlim([20,100]) 
        ylim([0.0,0.1])
        legend({'X','Y','T'})
    end
    
    
end

