
function [signal,fso] = gensigsources(type,par,N,path)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GENSIGSOURCES generates signals from sources 
% 
%type = {'senoide','wgn','audio',..., type_N}
%
%par = {par_1,par_2,...,par_N}:
%
%   type = 'senoide' -> par = [fso, duration, frequency]
%   type = 'audio' -> par = {path}
%   type = 'wgn' -> par = [fso, duration]
%
%N = Number of sources 
%   Note: the size of types = {...} and par = {...} 
%   must be N 
%
%path = path of the folder that has the signs  
%   Note: if filename is not empty, 
%   the function will generate the signals from 
%   the specified types and parameters 
%
%Examples:
%Generate 2 sinusoids with fs=24khz, 5s duration and 10hz frequency
%   GENSIGSOURCES({'senoide'},{[24000,5,10]},2,[]);
%
%Generate 10 audio signals read from 'path'
%   GENSIGSOURCES({'audio'},{'../../../_dataset/adult_female_speech.wav'},10,[]);
%
%Generate 2 mixed signals
%   GENSIGSOURCES({'senoide','audio'},{[24000,5,10],'../../../_dataset/adult_female_speech.wav'},2,[]);
%
%Generate signals read from the folder specified by 'path'   
%   GENSIGSOURCES([],[],[],['../../../_dataset/']);
%
%Generate 20 noise signals
%   GENSIGSOURCES({'wgn'},{[1e3,1]},[20],[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if ~isempty(path)          
        files  = dir([path '*.*']);
        signal = cell(1,length(files));
        fso = cell(1,N);
        if isempty(files)
            error(['No files in the folder: ' folder]) 
        end
        parfor i=3:length(files)
            file  = importdata([path  files(i).name]);                    
            if isstruct(file)
                if size(file.data) > 1
                    file.data = file.data(:,1);
                end
                signal{1,i-2} = file.data;
                fso{1,i-2} = file.fs;
            else
                signal{1,i-2} = file;
            end
        end
    else
        
        if (length(type) && length(par)) == 1
                type = repmat(type,1,N);
                par = repmat(par,1,N);
        end
        signal = cell(1,N);
        fso = cell(1,N);
%         for j=1:N
        parfor j=1:N
            switch type{j}
                case 'senoide'
                    duration = par{j}(2);
                    fso{1,j} = par{j}(1);
                    f = par{j}(3);
                    dt = 1/fso{1,j};
                    t = (0:dt:duration-dt)'; 
                    signal{1,j} = sin(2*pi*f*t);
                case 'wgn'
                    fso{1,j} = par{j}(1);   
                    duration = par{j}(2);
                    signal{1,j} = wgn(1,round((1+duration)*fso{1,j}),1); 
                case 'audio'
                    path = par{j};
                    file = importdata(path);
                    if isstruct(file)
                        if size(file.data) > 1
                            file.data = file.data(:,1);
                        end
                        signal{1,j} = file.data;
                        fso{1,j} = file.fs;
                    else
                        signal{1,j} = file;
                    end
                    
                otherwise
                    error(['Invalid signal type. Select an available type:  "senoide", "wgn", "audio".' ])
            end
        end
    
    end
end














