function  m_t = message_signal_generator(varargin)
% First input is signal type, second input is f_m, third input is time vector.
signal_type = varargin{1};
f_m = varargin{2};
t = varargin{3};
switch signal_type
    case 1
        m_t = (cos(2*pi*f_m*t));
    case 2
       m_t = (sawtooth(2*pi*f_m*t, 0.5));
    case 3
        m_t = (square(2*pi*f_m*t));
    case 4
        m_t = (cos(2*pi*f_m/2*t) + cos(2*pi*f_m*t));
end
end