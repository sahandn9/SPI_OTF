% Last updated on: 12/18/12

clc; clear; close all
%----------------------- Spatial domainm parameters -----------------------
x_size = 1; %Physical dimension
y_size = 1; %Physical dimension
N = 512;
dx = x_size/N;
dy = y_size/N;
x = 0:dx:1-dx;
y = 0:dy:1-dy;
[XX YY] = meshgrid(x ,y);
%--------------------- Frequency domain parameters ------------------------
fx_max = 1/(2*dx); % Nyquist frequency
fy_max = 1/(2*dy); % Nyquist frequency
dfx = fx_max/N;
dfy = fy_max/N;
fx = -fx_max/2:dfx:fx_max/2-dfx;
fy = -fy_max/2:dfy:fy_max/2-dfy;
[FXX FYY] = meshgrid(fx, fy);
%--------------------- Phase object parameters ----------------------------
% Sine
M =250; % Maximum number of cycles per x_size
lambda_x = 1:10:M*x_size; 
lambda_y = 1:10:M*y_size; 
psi_x = [0 2*pi./lambda_x]; % 0 is the DC term
psi_y = [0 2*pi./lambda_y];
w = 0.1; % Phase weakness constant

%==========================================================================

% ------------------------------- Phase Stepping --------------------------
phase_steps = [0 pi/4 pi/2 3*pi/4];
for p = 1:length(phase_steps)
    % ------------------ Spatial Phase Freq. Sweeping ---------------------
    for i=1:length(psi_y)
        for j=1:length(psi_x)
            
            Phi = cos(psi_x(j).*XX)+cos(psi_y(i).*YY);
            
            % Save the original Phi (no phase stepping needed)
            if (p==1)
               Phi_original{i,j} = Phi;
            end

            % Unity amplitude, phase(x,y) times the weakness factor times
            % the phase step
            phase_obj = 1.*exp(w*1i*Phi).*exp(1i.*phase_steps(p));

            % First lens (FFT)
            U1 = fftshift(fft2(phase_obj));

            % Fourier Plne (Filtering: Shearing in both X- and Y-direction)
            S = 1;
            shift_term_fx = exp(1i*2*pi.*FXX.*S/2) - exp(-1i*2*pi.*FXX.*S/2);
            shift_term_fy = exp(1i*2*pi.*FYY.*S/2) - exp(-1i*2*pi.*FYY.*S/2);

            Ud = (shift_term_fx + shift_term_fy).*U1;
            % Second lens (FFT)
            ud{p,i,j} = fft2(fftshift(Ud)); % 3-D data structure 
            % Detector (returns real values)
            img{p,i,j} = ud{p,i,j}.* conj(ud{p,i,j});      
        end
    end
    % ---------------------------------------------------------------------
 end
% -------------------------------------------------------------------------

% =========================================================================

% ------------------------- Phase Reconstruction --------------------------
for i=1:length(psi_y)
    for j=1:length(psi_x)
        beta = acos((img{1,i,j}-img{2,i,j}+img{3,i,j}-img{4,i,j})./...
            (2.*(img{2,i,j}-img{3,i,j})) );
        % RETURING COMPLEX VALUES!!!
        Phi_reconstructed{i,j} = atan(tan(beta/2).*...
            (img{1,i,j}-img{4,i,j}+img{2,i,j}-img{3,i,j})./...
            (img{2,i,j}+img{3,i,j}-img{1,i,j}-img{4,i,j}));        
    end
end
% -------------------------------------------------------------------------
figure;
imagesc(Phi_original{2,1}); colormap gray;
figure;
imagesc(Phi_reconstructed{2,1}); colormap gray;