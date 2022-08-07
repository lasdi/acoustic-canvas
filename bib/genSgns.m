function [Sgns]=genSgns(data_RIR, Fsbb,SNR, signal, fso,Nsig,name_dataset,pathDataset)

    USE_GPU = 0;
    if nargin < 7
        name_dataset =[];
        pathDataset =[];
    end
    if length(data_RIR) == 2
        RIRs= data_RIR{1,1};
        rot = data_RIR{1,2};
    else
        %Lê RIR da pasta (se data_RIR = {})
        folder = [pathDataset name_dataset '/RIRs/'];
        files  = dir([folder '*.mat']);
        if isempty(files)
            error(['No files in the folder: ' folder]) 
        end
        
        for j=1:length(files)
            data_rir= load([folder  files(j).name]);
            RIRs(:,j)= data_rir.h;
            rot(:,j) = data_rir.rot;
        end
    end
    
    [N,M] = size(RIRs);     
    
    %leitura do rotulo 
    [RT60, delay, dist_mic, fs] = labelread(rot,N);
    
    %Ajuste dos sinais das fontes
    [signal] = signaladjust(signal, fso, Fsbb, fs, N, Nsig);   
    
    %Núcleo de calculo dos sinais
    Sgns = zeros(Nsig,M);  
    parfor j = 1:M
%     for j = 1:M
        
        %Calcula o sinal que chega em cada mic j
        [rr] = calckernel(RIRs,signal,Nsig,USE_GPU, j);
        
        %Calcula a DRR (somente quando há 1 fonte)
        %%[DRR] = drr(signal, RIRs{1,j}, rr,N);
        DRR=[];

        %Insere ruido em cada mic j
        [Sgns(:,j)] = insertnoise(rr,SNR, fs, Fsbb);

        %Salva sinais de cada mic j (se name_dataset e pathdataset tiverem
        %o caminho)
        y = Sgns(:,j); 
        if ~ isempty(name_dataset)
        save_Sgn(y, SNR,RT60,delay,fs,dist_mic,DRR,j,N,name_dataset, pathDataset);        
        end
    end

end


function [RT60, delay, dist_mic, fs] = labelread(rot,N)
    if N < 2
        RT60=rot(1,1);
        delay = rot(2,2);
        dist_mic = rot(3,1); 
        fs = rot(4,1);
    else
        RT60=rot(1,1); delay = []; dist_mic = []; fs = rot(4,1);
    end
end


function [signal] = signaladjust(signal, fso, Fsbb, fs, N, Nsig)  
    parfor i=1:N
%     for i=1:N 
        aux = signal{1,i};        
        if fs ~= Fsbb
            aux = resample(aux,Fsbb,fso{1,i});
            aux = resample(aux,fs,Fsbb);
        else
            aux = resample(aux,fs,fso{1,i});
        end
        aux = aux - mean(aux);
        aux = aux/max(abs(aux));
        aux = aux(1:Nsig);
        signal{1,i} = aux;
    end
end


function [rr] = calckernel(RIRs,signal,Nsig,USE_GPU,j)

    [N,~] = size(RIRs);
    rr = zeros(N,Nsig);
    
    if (USE_GPU==1)
%         for i = 1:N
        parfor i = 1:N
            rir = RIRs{i,j};
            A_g = gpuArray(signal{1,i});
            B_g = gpuArray(rir);
            rr_g = conv(A_g, B_g);
            r_conv = gather(rr_g);
            rr(i,:) = r_conv(1:Nsig);
        end
    else
%         for i = 1:N
        parfor i = 1:N
            rir = RIRs{i,j};
            A_g = signal{1,i};
            B_g = rir(:);
            r_conv = conv(A_g, B_g);
            rr(i,:) = r_conv(1:Nsig);   
        end
    end
    rr = sum(rr,1);    
end



function [DRR] = drr(signal, rir, rr,N)
    if N < 2
        coef_att = max(rir);
        Sgn_d = cell2mat(signal)*coef_att;
        Pd = mean(Sgn_d.^2);
        Pt = mean(rr.^2);
        DRR = pow2db(Pd/(Pt-Pd));
    else
        DRR=[];
    end

end


function [Sgns] = insertnoise(rr,SNR, fs, Fsbb)
    
    Nsig = length(rr);
    OSR = fs/Fsbb;
    pwr_sin = pow2db(bandpower(rr));
    pwr_noise = pwr_sin-SNR + 10*log10(OSR);
    
    w = wgn(Nsig,1,pwr_noise,'dBW');
    Hd = lpf_iir_gen(fs);
    
    if OSR~=1
        w = filter(Hd, w);
    end
    Sgns = rr(:)+ w(:);
end


function save_Sgn(y,SNR,RT60,delay,fs,dist_mic,DRR,i,N,name_dataset, path_dataset)

    if ~isempty(name_dataset)
        name_Folder = name_dataset;        
        path = fullfile(path_dataset,name_Folder);
        filename = ['Sgn_' num2str(SNR) '_' num2str(RT60*1000) '_' num2str(fs) '_' num2str(i-1) '.mat'];
        if ~exist(path)
            mkdir(path);
        end
        save([path '/' filename],'y')
        
        if N < 2
            fid = fopen([path '/' filename(1:end-6) '.rot'],'w');
            rot = [RT60;SNR;delay; dist_mic;fs;DRR(:)];
            fprintf(fid,'%d\n',rot);
            fclose(fid);
        end
    end
end 








