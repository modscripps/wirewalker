function dnum = unixtime2datenum(utime)

% unixtime2datenum - Convert unix time (i.e. the number of seconds
%                    since 1970-01-01) MATLAB datenum format
%
% Syntax:  dnum = unixtime2datenum(utime)
% 
% Converts "unix time", (i.e. POSIX time, the number of seconds since
% 1970-01-01 00:00:00.000) to MATLAB datenum format.
%
% Inputs:
%    utime - Unix time. The number of seconds since 1970-01-01.
% 
% Outputs:
%    dnum - MATLAB datenum
%
% Example: 
%    datestr(unixtime2datenum(1.420070400000000e+09))
%
%    ans =
%
%         01-Jan-2015
%
% See also: datenum2unixtime, RSKtime2datenum, datenum2RSKtime
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com

dnum = datenum(1970,1,1,0,0,utime)';
