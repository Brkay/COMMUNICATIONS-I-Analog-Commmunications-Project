function [M_f, x_t, X_f, S_f, B_exp] = fm_generator(varargin)
% First variable is m_t.
% Second variable is k_f.
% Third variable is T value.
% Fourth variable is t_s.
% Fifth variable is f_c.
% Sixth variable is A_c.
% Seventh variable is N_f.


%For the purpose of project, we gave just first two inputs in our main script
%and made rest of the variables  same with given values in the project,
% but they can be changed via giving different inputs. In order to have a
% correct analysis, T and t_s values should be given accordingly. (They
% should be same when m_t was constructed, so we took it as it was stated
% in main script; however, for different size of inputs(time vectors or
% sampling rates) we still leave it as variable.
m_t = varargin{1};
k_f = varargin{2};
% These variables have default  values, but if input is given different we
% change it in switch statement.
T = 1;
f_c=20e3;
t_s=1e-6;
N_f = 5e6;
A_c = 1;
switch nargin
    case 3
        T = varargin{3};
    case 4
        T = varargin{3};
        t_s = varargin{4};
    case 5
        T = varargin{3};
        t_s = varargin{4};
        f_c = varargin{5};
        
    case 6
        T = varargin{3};
        t_s = varargin{4};
        f_c = varargin{5};
        A_c = varargin{6};
    case 7
        T = varargin{3};
        t_s = varargin{4};
        f_c = varargin{5};
        A_c = varargin{6};
        N_f = varargin{7};
end
N = T / t_s;
f_s=(N/N_f)*(1/T);
t = (0:t_s:T-t_s).';
F = 1/t_s;
f=((-F/2):f_s:((F/2)-f_s)).';

length_of_message_signal = length(m_t); % Length of the message signal.

m_t_integral = cumsum(m_t)*t_s; % Integrated signal.

M_f= fft(m_t,N_f)./ length_of_message_signal;
M_f= fftshift(M_f); % Two sided fft of message signal.

x_t=A_c* cos(2*pi*f_c*t+2*pi*k_f*m_t_integral); % FM modulated signal.

X_f= fft(x_t,N_f)./ length_of_message_signal;
X_f= fftshift(X_f); % Two sided fft of FM modulated signal.

S_f=(abs(X_f).^2)/T; % PSD of FM modulated signal.

P=sum(S_f)*f_s/2; % One side of the spectrum is examined.
B_exp=0; % Initial guess for Bexp is 0.
P_eff=0;
% Increase symmetrically indexes to achieve 98% power
% and find effective bandwidth.
S_b=S_f;
[~, power_frequency_fc] = min ( abs( f - f_c ) );
while (P_eff/P) < 0.98 % Power is symmetrical to fc frequency.
    s_b_1= power_frequency_fc - B_exp/ (2*f_s);
    s_b_2 = power_frequency_fc + B_exp/ (2*f_s);
    P_eff=sum(S_b(s_b_1:s_b_2))*f_s;
    B_exp=B_exp+ 1;
end
end