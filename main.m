C = [0.5 1 2 4];
N = 6;
tt = 10;
Fs = 1000;
T = 2;
no_functions = 5;

A = zeros(N);
B = zeros(N);

%=================================PART 1=================================
%================== Generation and Plotting of PSWFs ====================
%   use command "lambda" to display the eigen values of the PSWFs
for c_idx=1:4
    c = C(c_idx);
    OMEGA = 2*c/T;
    
    t = -tt:1/Fs:tt;
    mid_t_idx = floor(Fs*(tt-T/2))+1:floor(Fs*(tt+T/2));

    [psi, lambda] = pswf(c,T,tt,Fs,N);

    A(c_idx,:) = lambda;
    
    figure;
    hold on;
    no_graph=6;
    for n=0:2:no_graph-1
        plot(t,real(psi(n+1,:)),'Color',[n/no_graph, n/no_graph, 1], 'LineWidth', 2);
    end
    xlabel("x");
    ylabel("\psi_n(x)");
    saveas(gcf, "report/c"+ c + "_e.jpg");
    hold off;

    figure;
    hold on;
    no_graph=6;
    for n=1:2:no_graph-1
        plot(t,real(psi(n+1,:)),'Color',[1, n/no_graph,n/no_graph], 'LineWidth', 2);
    end
    xlabel("x");
    ylabel("\psi_n(x)");
    hold off;
    saveas(gcf, "report/c"+ c + "_o.jpg");
    
end

%========================= Part 2 =======================================
%============ Determining inner product of all PSWFs ====================
%   use command "A" for inner product over finitely long(inf) time
%   use command "B" for inner product over finite time

A = real(psi(1:no_functions,:) * psi(1:no_functions,:)' / Fs);
B = real(psi(1:no_functions,mid_t_idx) * psi(1:no_functions,mid_t_idx)' / Fs);

%========================= Part 3 =======================================
%============ applying BD operation and plotting ========================
BD_psi = zeros(size(psi));

for n=1:size(t,2)
    x = t(n);
    sinc = sin(OMEGA*(x-t(mid_t_idx))) ./( pi*(x-t(mid_t_idx)));
    sinc(isnan(sinc)) = OMEGA/pi;
    BD_psi(:,n) = psi(:,mid_t_idx) * sinc' ./Fs;
end

for n=1:no_functions
    figure;
    hold on;
    plot(t,BD_psi(n,:),  'LineWidth', 2);
    plot(t,lambda(1,n)*psi(n,:),  'LineWidth', 2);
    xlabel("x");
    legend("BD \psi_"+ string(n-1) + "(x)", "\lambda_" + string(n-1) + "\psi_" + string(n-1) + "(x)");
    saveas(gcf, "report/BD_"+ n + ".jpg");
    hold off;
end


%=======================Part 4===========================================
%============== PSWF basis expansion of L2 sinc signal ==================
%==================and corresponding plots ==============================
%   use command "energy_ratio_tl", "energy_ratio" to display the energy 
%   ratios in finite and infinite domain
signal_t_lim = zeros(1, size(t,2));
signal_t_lim(mid_t_idx) = sin(pi*t(mid_t_idx)) ./ (pi* t(mid_t_idx));
signal_t_lim( isnan(signal_t_lim)) = 1;

signal = sin(pi*t) ./ (pi* t);
signal( isnan(signal)) = 1;
    
energy_ratio_tl = zeros(1, size(C,2));
energy_ratio = zeros(1, size(C,2));

gamma_n_tl = zeros(size(C,2), no_functions);
gamma_n = zeros(size(C,2), no_functions);

for c_idx=1:4
    c = C(c_idx);
    OMEGA = 2*c/T;
    
    t = -tt:1/Fs:tt;
    mid_t_idx = floor(Fs*(tt-T/2))+1:floor(Fs*(tt+T/2));

    [psi, lambda] = pswf(c,T,tt,1000,N);
    
    gamma_n_tl(c_idx,:) = signal_t_lim * psi(1:no_functions,:)' ./Fs;
    gamma_n(c_idx,:) = signal * psi(1:no_functions,:)' ./Fs;

    reconstructed_t_lim_signal = gamma_n_tl(c_idx,:) * psi(1:no_functions,mid_t_idx);
    reconstructed_infinite_signal = gamma_n(c_idx,:) * psi(1:no_functions,:);

    energy_ratio_tl(1, c_idx) = (rms(reconstructed_t_lim_signal) / rms(signal_t_lim))^2;
    energy_ratio(1, c_idx) = (rms(reconstructed_infinite_signal) / rms(signal))^2;
    
    figure;
    hold on;
    plot(t(mid_t_idx), signal_t_lim(mid_t_idx), 'LineWidth', 2);
    plot(t(mid_t_idx),reconstructed_t_lim_signal, 'LineWidth', 2);
    legend("Original", "Reconstructed");
    xlabel("time");
    axis([-T/2, T/2, 0, 1]);
    saveas(gcf, "report/recon_tl_"+ c + ".jpg");
    hold off;
    
    figure;
    hold on;
    plot(t, signal, 'LineWidth', 2);
    plot(t,reconstructed_infinite_signal, 'LineWidth', 2);
    xlabel("time");
    legend("Original", "Reconstructed");
    saveas(gcf, "report/recon_"+ c + ".jpg");
    hold off;
    
end

