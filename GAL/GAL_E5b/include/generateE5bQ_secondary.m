function E5bQ_secondary = generateE5bQ_secondary(PRN)
% This function generates Galileo E1B code in bipolar format (-1, +1)
% PRNs 1 to 50 (no care is taken if PRN number is > 50)
% The memory codes are stored in the file gale1bcode.dat
%
% E1Bcode = generateE1Bcode(PRN)
%
%   Inputs:
%       PRN             - PRN number of the sequence.
%
%   Outputs:
%       E5bQ_secondary  - a vector containing the desired E5BQ secondary
%                          code sequence(chips).
%
%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: generateE5BIcode.m,v 1.1.2.5 2017/04/24 22:00:00 dpl Exp $

%% Secondary code table from Galileo-OS-SIS-ICD (2016)
secondary_code = [...
    'CFF914EE3C6126A49FD5E5C94',    'FC317C9A9BF8C6038B5CADAB3',...
    'A2EAD74B6F9866E414393F239',    '72F2B1180FA6B802CB84DF997',...
    '13E3AE93BC52391D09E84A982',    '77C04202B91B22C6D3469768E',...
    'FEBC592DD7C69AB103D0BB29C',    '0B494077E7C66FB6C51942A77',...
    'DD0E321837A3D52169B7B577C',    '43DEA90EA6C483E7990C3223F',...
    '0366AB33F0167B6FA979DAE18',    '99CCBBFAB1242CBE31E1BD52D',...
    'A3466923CEFDF451EC0FCED22',    '1A5271F22A6F9A8D76E79B7F0',...
    '3204A6BB91B49D1A2D3857960',    '32F83ADD43B599CBFB8628E5B',...
    '3871FB0D89DB77553EB613CC1',    '6A3CBDFF2D64D17E02773C645',...
    '2BCD09889A1D7FC219F2EDE3B',   '3E49467F4D4280B9942CD6F8C',...
    '658E336DCFD9809F86D54A501',    'ED4284F345170CF77268C8584',...
    '29ECCE910D832CAF15E3DF5D1',    '456CCF7FE9353D50E87A708FA',...
    'FB757CC9E18CBC02BF1B84B9A',    '5686229A8D98224BC426BC7FC',...
    '700A2D325EA14C4B7B7AA8338',    '1210A330B4D3B507D854CBA3F',...
    '438EE410BD2F7DBCDD85565BA',    '4B9764CC455AE1F61F7DA432B',...
    'BF1F45FDDA3594ACF3C4CC806',    'DA425440FE8F6E2C11B8EC1A4',...
    'EE2C8057A7C16999AFA33FED1',    '2C8BD7D8395C61DFA96243491',...
    '391E4BB6BC43E98150CDDCADA',    '399F72A9EADB42C90C3ECF7F0',...
    '93031FDEA588F88E83951270C',    'BA8061462D873705E95D5CB37',...
    'D24188F88544EB121E963FD34',    'D5F6A8BB081D8F383825A4DCA',...
    '0FA4A205F0D76088D08EAF267',    '272E909FAEBC65215E263E258',...
    '3370F35A674922828465FC816',    '54EF96116D4A0C8DB0E07101F',...
    'DE347C7B27FADC48EF1826A2B',    '01B16ECA6FC343AE08C5B8944',...
    '1854DB743500EE94D8FC768ED',    '28E40C684C87370CD0597FAB4',...
    '5E42C19717093353BCAAF4033',    '64310BAD8EB5B36E38646AF01'];

% To prevent data overflow, use two segemnts to perform converting
% hexadecimal to decimal
first_segement = zeros(1,13*4);
second_segement = zeros(1,12*4);

% Take hexadecimal secondary code for # PRN
code = secondary_code(25*(PRN-1)+1:25*PRN);

% Perform convertion from hexadecimal to binary for segment 1
code_temp = dec2bin(base2dec(code(1:13),16))-48;
first_segement(end-length(code_temp)+1:end) = code_temp;

% Perform convertion from hexadecimal to binary for segment 2
code_temp = dec2bin(base2dec(code(14:end),16))-48;
second_segement(end-length(code_temp)+1:end) = code_temp;

% Joint the two segments
E5bQ_secondary = [first_segement second_segement];

% Convert to bipolar format (0-->1; 1-->-1)
E5bQ_secondary = 1 - 2* E5bQ_secondary;

