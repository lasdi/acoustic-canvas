function S = steeredResponseCBF(inputSignal, fs, elementWeights, xPos, yPos,...
    zPos, f, c, thetaScanAngles, phiScanAngles,NFFT,f_min,f_max,...
    scanningPointsX,scanningPointsY,scanningPointsZ)
%STEEREDRESPONSECBF Summary of this function goes here
%   Detailed explanation goes here

N = length(xPos);
% L = length(thetaScanAngles)*length(phiScanAngles);
L = numel(thetaScanAngles);
theta = reshape(thetaScanAngles,1,L);
phi = reshape(phiScanAngles,1,L);

f_min = floor(f_min*NFFT/fs);
if f_min==0
    f_min=1;
end
f_max = ceil(f_max*NFFT/fs);
f_v = ((f_min:f_max)-1) * fs/NFFT;

%% Cross Spectral Matrix
R = zeros(N,N,NFFT/2+1);
for i = 1:N
    parfor j = 1:N
        if i==j
            continue;
        end
        R(i,j,:) = cpsd(inputSignal(i,:),inputSignal(j,:),NFFT);
%         R(j,i,:) = -R(i,j,:); 
    end
end

%% Steering Vector
% tic;disp("SV")
scanningPointsX=scanningPointsX(:);
scanningPointsY=scanningPointsY(:);
scanningPointsZ=scanningPointsZ(:);

% % Sarradj Formulation I
% h = zeros(N,L,NFFT);
% parfor l = 1:L
%     d = sqrt( (scanningPointsX(l)-xPos).^2 + ...
%                  (scanningPointsY(l)-yPos).^2 + ...
%                  (scanningPointsZ(l)-zPos).^2 );
%     d=d-d(1);
%     k = 2*pi*f_v/c;
%     for ff = f_min:f_max
%         h(:,l,ff) = exp(-1j*k(ff-f_min+1)*d)/N;
%     end
% end


%% Teste
h = zeros(N,L,NFFT);
Z = zeros(L,NFFT);
for l = 1:L
    d = sqrt( (scanningPointsX(l)-xPos).^2 + ...
                 (scanningPointsY(l)-yPos).^2 + ...
                 (scanningPointsZ(l)-zPos).^2 );
    d_r = d(1)./d;
    d = d-d(1);
    k = 2*pi*f_v/c;
    for ff = f_min:f_max
        k_aux = k;
        g = d_r.*exp(-1j*k_aux(ff-f_min+1)*d);
        h_temp = (g./(sqrt(N)*sqrt(g*g'))).';
        Z(l,ff)= h_temp'*R(:,:,ff)*(h_temp);
    end
end
%%


% % Sarradj Formulation IV
% h = zeros(N,L,NFFT);
% for l = 1:L
%     d = sqrt( (scanningPointsX(l)-xPos).^2 + ...
%                  (scanningPointsY(l)-yPos).^2 + ...
%                  (scanningPointsZ(l)-zPos).^2 );
%     d_r = d(1)./d;
%     d = d-d(1);
%     k = 2*pi*f_v/c;
%     for ff = f_min:f_max
%         g = d_r.*exp(-1j*k(ff-f_min+1)*d);
%         h(:,l,ff) = g./(sqrt(N)*sqrt(g*g'));
%     end
% end
% 
% % toc
% %% CBF
% Z = zeros(L,NFFT);
% 
% for l = 1:L
%     for ff = f_min:f_max
%         h_temp = h(:,l,ff);
%         Z(l,ff)= h_temp'*R(:,:,ff)*(h_temp);
%     end
% end

S = sum((abs(Z(:,f_min:f_max))),2);
% S = S/max(S);
% S = 10*log10(S);
S = reshape(S,size(thetaScanAngles));

end

