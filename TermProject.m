%% EE-435 Term Project, Berkay Yaldiz 2232940, Melih Can Zerin 2233088
%% Part 1
% Part b
% We found at part a) that we need to sum all the rectangle heights and
% multiply it with interval delta Tau.
% Let x(t) = t^2 and T = 10.
clear
close all
T = 10;
a = 0;
b = T;
n = 1000000;% As n increases error gets smaller.
delta = (b-a)/n;
t = a:delta:b;
X = t.*t;
Area = cumsum(X) *delta ;
p=cumtrapz(t,X); % Matlab's built in function.
error=Area-p;
clearvars -except error
% Part c, vi
% Since we just need a 2pi interval for DFT frequencies, after taking fft; fftshift can be
% used to represent DFT outputs between -pi and pi intervals (Sampling frequency normalized to 2*pi).
%% Part 2
% We wrote according to given conventions. We used column vectors for
% signals.
% message_signal_generator.m and fm_generator.m functions are the answers for this part.
Ac = 1;
T = 1;
t_s = 1e-6;

N = T / t_s;

N_f = 5e6;

f_s=(N/N_f)*(1/T);
F = 1/t_s;

f=((-F/2):f_s:((F/2)-f_s)).';

t = (0:t_s:T-t_s).';

f_c = 20e3;

%% Part 3
signal_number = 7; % Number of signals can be changed here, but f_m and k_f inputs should be entered by hand!
signal_type_part3 = 1; % Signal type input ( 1 for cos, 2 for sawtooth, 3 for square, 4 for sum of 2 cosine signals)
f_m = [ repmat(1e3,1,3) repmat(2e3, 1, 2) repmat(1e3,1,2) ]; % 7 different signal parameters in part 3 (a-b-c-d-e-f-g)
k_f  = [ 0.1e3 1e3 10e3 1e3 2e3 10e3 10e3];
theoretical_beta = zeros(1, signal_number); % Theoritical beta values.
B_c = zeros(1,signal_number);
% Now we will construct signals with our function in a for loop, each
% column represents a, b, c, d, e, f, g signals respectively.

% Initialize used variables for speed.
length_of_signal = length(t);
m_t = zeros(length_of_signal ,signal_number);
x_t = zeros(length_of_signal , signal_number);

M_f = zeros(N_f, signal_number);
X_f = zeros(N_f, signal_number);
S_f = zeros( N_f, signal_number);
P = zeros(1, signal_number);
P_eff= zeros(1, signal_number);
B_exp=zeros(1, signal_number);
for ii=1: signal_number
    if isequal(ii,6)
        signal_type_part3 = 2; % Make message signal sawtooth.
    elseif isequal(ii,7)
        signal_type_part3 = 3; % Make message signal square.
    end
    m_t(:,ii) = message_signal_generator(signal_type_part3,f_m(ii),t); % First input is signal type, second input is f_m, third input is time vector.
    [M_f(:,ii), x_t(:,ii), X_f(:,ii), S_f(:,ii), B_exp(:,ii)] = fm_generator(m_t(:,ii), k_f(ii) );
    
    theoretical_beta(ii) = k_f(ii) / f_m(ii); % A_m is 1.
    B_c(ii) = 2* f_m(ii) * (theoretical_beta(ii)  + 1);
    
    % Take interval between  -2*fm and 2*fm for plotting M(f) and X(f) and
    % 0 - 2/f_m for time plot index.
    [~, index_time] = min( abs ( t - ( 2 / f_m(ii) ) ) );
    if  ii < 6
        [~, M_f_index1] = min( abs ( f - (-2*f_m(ii) ) ) ); % We used min-abs approach because there was a rounding-precision problem while determining the index.
        [~, M_f_index2] = min( abs ( f - (2*f_m(ii) ) ) ) ;
    else
        [~, M_f_index1] = min ( abs( f - (-5 * f_m(ii) ) ) );
        [~, M_f_index2] = min ( abs ( f - (5*f_m(ii) ) ) );
    end
    
    [~, X_f_index1] = min ( abs( f - ( f_c - max( B_exp(ii), B_c(ii) ) ) ) );
    [~, X_f_index2] = min ( abs ( f - ( f_c + max( B_exp(ii), B_c(ii) ) ) ) );

    figure('Position',[0 0 1920 1080]) % In order to plot figures with full size.
    subplot(2,2,1)
    plot(t(1:index_time), m_t(1:index_time,ii) )
    grid on
    if ii < 6
        title({'Time plot of message signal (Cosinus)', ['frequency = ',num2str(f_m(ii)) ,' Hz'] } )
    elseif isequal(ii,6)
        title({'Time plot of message signal (Sawtooth)', ['frequency = ',num2str(f_m(ii)) ,' Hz'] } )
    else
        title({'Time plot of message signal (Square)', ['frequency = ',num2str(f_m(ii)) ,' Hz'] } )
    end
    ylabel('Amplitude (V)')
    xlabel(' Time (seconds)')
    
    subplot(2,2,2)
    plot(t(1:index_time), x_t(1:index_time,ii) )
    grid on
    title({ 'Time plot of modulated signal', [ 'f_m =',num2str(f_m(ii)) ,' Hz', ' k_f = ',num2str(k_f(ii)), ' Hz/V' ] } )
    ylabel('Amplitude (V)')
    xlabel(' Time (seconds)')
    
    subplot(2,2,3)
    plot( f ( M_f_index1:M_f_index2), abs( M_f ( M_f_index1:M_f_index2,ii) ) )
    grid on
    if ii < 6
        title( { 'Frequency plot of message signal' , ['f_m = ',num2str(f_m(ii)) ,' Hz'] } )
    elseif isequal(ii,6)
        title( { 'Frequency plot of message signal' , ['f_m = ',num2str(f_m(ii)) ,' Hz'] } )
    else
        title( { 'Frequency plot of message signal' , ['f_m = ',num2str(f_m(ii)) ,' Hz'] } )
    end
    ylabel('Amplitude (V)')
    xlabel(' Frequency (Hz)')

    subplot(2,2,4)
    plot( f ( X_f_index1:X_f_index2), abs(X_f ( X_f_index1:X_f_index2) ) ) ;
    grid on
    title( { 'Frequency plot of modulated signal',[' f_m = ',num2str(f_m(ii)) ,' Hz k_f = ', num2str(k_f(ii)), ' Hz/V'] } )
    ylabel('Amplitude (V)')
    xlabel(' Frequency (Hz)')
end
%% Part 4
k_f4 = 10e3;
f_m4 = 100;
signal_number_part4 = 3;

m_t_part4 = zeros(length_of_signal, signal_number_part4);
B_c_part4 = zeros(1,signal_number_part4);
theoretical_beta_part4 = zeros(1,signal_number_part4);
S_f_part4 = zeros(N_f, signal_number_part4);
B_exp_part4 = zeros(1,signal_number_part4);
for ii=1 : signal_number_part4
    m_t_part4(:,ii) = message_signal_generator(ii, f_m4, t); % ii=1 - cos,  ii=2 -sawtooth, ii=3 square. If signal number changes, please change the input to message signal genetor accordingly.
    [~,~,~,S_f_part4(:,ii), B_exp_part4(ii)] = fm_generator(m_t_part4(:,ii), k_f4);% We do not need x_t, X_f, M_f for this part.
    
    theoretical_beta_part4(ii) = k_f4 / f_m4; % A_m is 1.
    B_c_part4(ii) = 2* f_m4 * (theoretical_beta_part4(ii)  + 1);
    
    figure('Position',[0 0 1920 1080])
    subplot(3,1,1)
    [~,time_index] =min ( abs ( t - 2 ./ f_m4) );
    plot( t (1 : time_index), m_t_part4 (1: time_index,ii) ) % Plot between 0 and 2/f_m
    if isequal(ii,1)
        title( { 'Time plot of message signal (Cosinus)', ['frequency =',num2str(f_m4) ,' Hz'] } )
    elseif isequal(ii,2)
        title( { 'Time plot of message signal (Sawtooth)', ['frequency =',num2str(f_m4) ,' Hz'] } )
    else
        title( { 'Time plot of message signal (Square)', ['frequency =',num2str(f_m4) ,' Hz'] } )
    end
    ylabel('Amplitude (V)')
    xlabel(' Time (seconds)')
    grid on
    subplot(3,1,2)
    [~, freq_index_part4_1] = min ( abs( f - ( f_c - max( B_exp_part4(ii), B_c_part4(ii) ) ) ) );
    [~, freq_index_part4_2] = min ( abs( f - ( f_c + max( B_exp_part4(ii), B_c_part4(ii) ) ) ) );
    plot( f (freq_index_part4_1 : freq_index_part4_2 ),S_f_part4(freq_index_part4_1 : freq_index_part4_2,ii)  )
    if isequal(ii,1)
        title( { 'PSD of modulated signal (Cosinus)', ['frequency =',num2str(f_m4) ,' Hz'] } )
    elseif isequal(ii,2)
        title( { 'PSD plot of modulated signal (Sawtooth)', ['frequency =',num2str(f_m4) ,' Hz'] } )
    else
        title( { 'PSD plot of modulated signal (Square)', ['frequency =',num2str(f_m4) ,' Hz'] } )
    end
    ylabel('Amplitude (V)')
    xlabel(' Frequency (Hz)')
    grid on
    subplot(3,1,3)
    myhistogram=histogram( k_f4 * m_t_part4(:,ii));
    myhistogram.Normalization='pdf'; % Since histogram is generally used to calculate pdf of a signal, we normalized it.
    grid on
    if isequal(ii,1)
        title( { 'Histogram plot of message signal (Cosinus)', ['frequency =',num2str(f_m4) ,' Hz'] } )
    elseif isequal(ii,2)
        title( { 'Histogram plot of message signal (Sawtooth)', ['frequency =',num2str(f_m4) ,' Hz'] } )
    else
        title( { 'Histogram plot of message signal (Square)', ['frequency =',num2str(f_m4) ,' Hz'] } )
    end
    
end
%% Part 5
% Part a
% We can write FM signal expression for NBFM as follows ( EE435_angle
% modulation)
% x_t_NBFM_c=cos(2*pi*f_c*t)-((sin(2*pi*f_c*t))*2*pi*k_f_5_NBFM.*m_t_5_integral);

% Part b
f_m_5=1e3;
m_t_5 = message_signal_generator(1,f_m_5,t);
beta_5_b=0.2;
k_f_5_b=beta_5_b*f_m_5;
[~,x_t_5_b,~,~,~]=fm_generator(m_t_5,k_f_5_b);

beta_5_NBFM=0.1; % Since NBFM beta is half of the required beta, we selected NBFM_fc as fc/2 (10kHz).
f_c_NBFM1 = f_c/2;
k_f_5_NBFM=beta_5_NBFM*f_m_5;
x_t_NBFM_b=cos(2*pi*f_c_NBFM1*t)-beta_5_NBFM*(sin(2*pi*f_c_NBFM1*t)) .* sin(2*pi*f_m_5*t);
x_armstrong_t_b=2*((x_t_NBFM_b.^2)-0.5);
[~,time_index_part5] = min( abs ( t - ( 2 / f_m_5 ) ) );
figure
plot( t (1:time_index_part5), x_t_5_b( 1:time_index_part5) )
hold on;
plot(t (1:time_index_part5) , x_armstrong_t_b(1:time_index_part5));
title('Direct Generated and Armstrong''s signals (beta =0.2) ')
ylabel('Amplitude (V)')
xlabel('Time (seconds)')
legend('Direct Generated','Armstrong')
hold off;
grid on
error_for_NBFM1 = mean((abs(x_t_5_b-x_armstrong_t_b)));
% Mean error for part b can be seen above (beta = 0.2).

% Part c
beta_5_c=12.8;
k_f_5_c=beta_5_c*f_m_5;
[~,x_t_5_c,~,~,~]=fm_generator(m_t_5,k_f_5_c);

beta_5_NBFM2 = 0.1; 
f_c_NBFM2 = f_c/(2^7); % % Since NBFM beta is 1/128 of the required beta, we selected NBFM_fc as fc/128 (0.1563) kHz.
k_f_5_NBFM2=beta_5_NBFM2*f_m_5;
x_t_NBFM_c=cos(2*pi*f_c_NBFM2*t)-beta_5_NBFM2*(sin(2*pi*f_c_NBFM2*t)) .* sin(2*pi*f_m_5*t);
x_armstrong_t_c = x_t_NBFM_c;
for ii = 1:7
    x_armstrong_t_c=2*((x_armstrong_t_c.^2)-0.5);
end
figure
plot( t(1:time_index_part5), x_t_5_c(1:time_index_part5) ) % Plot from 0 to 2/fm.
hold on;
plot( t(1:time_index_part5) , x_armstrong_t_c(1:time_index_part5) ) 
hold off;
title('Direct Generated and Armstrong''s signals (beta =12.8) ')
ylabel('Amplitude (V)')
xlabel('Time (seconds)')
hold off;
grid on
legend('Direct Generated','Armstrong')
error_for_NBFM2 = mean((abs((x_t_5_c-x_armstrong_t_c))));
% Mean error for part b can be seen above (beta = 12.8).
% As we can see from the plot, Armstrong's method catches the original
% signal after nearly 1ms.
%% Part 6
% Part a
f1_part6 = 18e3;
f2_part6 = 22e3;
Ts = 250e-6;
A_c_6 = 1;
f_c_6 = 20e3;
kf_6 = 2e3;
t_6 = (0:t_s:(Ts-t_s)).';
N_s = 4e3;
realization_number = 1000;
X_f_part6 = zeros( N_f,1);
symbol_length=length(t_6);
x_t_6 = zeros( N_s*length(t_6),1);
m_t_6 = zeros( N_s*length(t_6), 1);
X_f_6 = zeros( N_f, 1); % FFT length is N_f = 5e6.
S_f_6_temp = zeros( N_f, 1);
a_k = randi( [0 1], N_s, realization_number);
zero_indexes = a_k ==0;
a_k(zero_indexes) = -1;
time_vector_part6 = zeros(N_s*length(t_6),1);
% Consider the following lines for speed.
for kk = 1: N_s
time_vector_part6( symbol_length*(kk-1) + 1 : (kk) *symbol_length) = ( kk-1).*Ts + t_6; % construct time vector for each a_k.
end % Construct time vector just once.

symbol_1 = (f1_part6 - f_c_6) ./ kf_6 *ones(symbol_length,1);
symbol_2 = (f2_part6 - f_c_6) ./ kf_6 *ones(symbol_length,1); 
length_of_signal_part6 = length(x_t_6);
% Calculate just one time these three parameters.

for ii = 1:realization_number
    for kk = 1: N_s
        if a_k(kk,ii) == -1
            x_t_6(symbol_length*(kk-1) + 1 : (kk) *symbol_length) =  cos( 2 * pi * f1_part6 * time_vector_part6(symbol_length*(kk-1) + 1 : (kk) *symbol_length) );
            if ii == 1000 % Consider only one realization of message signal. (Last one, also true for modulated signal.)
                m_t_6(symbol_length*(kk-1) + 1 : (kk) *symbol_length) = symbol_1 ;
            end
        else
            x_t_6(symbol_length*(kk-1) + 1 : (kk) *symbol_length) =   cos( 2 * pi * f2_part6 * time_vector_part6(symbol_length*(kk-1) + 1 : (kk) *symbol_length) );
            if ii == 1000 
                m_t_6(symbol_length*(kk-1) + 1 : (kk) *symbol_length) = symbol_2;
            end
        end
    end
    X_f_part6 = fft(x_t_6,N_f) ./ length_of_signal_part6 ;  % Since trials will be long, we did not use our function and directly calculated N_f points FFT.
    S_f_6_temp = S_f_6_temp + (abs(X_f_part6).^2); 
end
time_index_part6 = 1.5e-3 / t_s;
figure
subplot(3,1,1)
plot( t(1:time_index_part6), m_t_6(1:time_index_part6) )
title('Message signal for FSK')
xlabel('Times (Seconds)')
ylabel('Amplitude (V)')
grid on

subplot(3,1,2)
plot( t(1:time_index_part6), x_t_6(1:time_index_part6) )
title('FSK modulated signal')
xlabel('Times (Seconds)')
ylabel('Amplitude (V)')
grid on

[~, freq_index_part6_1] = min ( abs( f - 8e3 ) );
[~, freq_index_part6_2] = min ( abs( f -  32e3  ) );
S_f_6_temp = fftshift(S_f_6_temp);
S_f_6 = S_f_6_temp / T / realization_number;
subplot(3,1,3)
plot( f ( freq_index_part6_1: freq_index_part6_2), 10*log10(S_f_6(freq_index_part6_1:freq_index_part6_2)))
title('PSD of FSK modulated signal in dB scale')
xlabel('Frequency (Hz)')
ylabel('Amplitude (V in dB)')
grid on
 
% Part iii- Calculation of Bexp:
P_6 = sum(S_f_6 )*f_s /2; % One sided Power.
P_eff_6 = 0;
B_exp_6 = 0;
[~, power_frequency_fc] = min ( abs( f - f_c ) );
S_b=S_f_6;
while (P_eff_6/P_6)<0.98 % Power is symmetrical to fc frequency.
    s_b_part6_1= power_frequency_fc - B_exp_6/(2*f_s); 
    s_b_part6_2 = power_frequency_fc + B_exp_6 / (2*f_s);
    P_eff_6=sum(S_b(s_b_part6_1:s_b_part6_2))*f_s;
    B_exp_6=B_exp_6+1;
end






