function E5aQ_secondary = generateE5aQ_secondary(PRN)
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
secondary_code = [ ...
    '83F6F69D8F6E15411FB8C9B1C',    '66558BD3CE0C7792E83350525', ...
    '59A025A9C1AF0651B779A8381',    'D3A32640782F7B18E4DF754B7', ...
    'B91FCAD7760C218FA59348A93',    'BAC77E933A779140F094FBF98', ...
    '537785DE280927C6B58BA6776',    'EFCAB4B65F38531ECA22257E2', ...
    '79F8CAE838475EA5584BEFC9B',    'CA5170FEA3A810EC606B66494', ...
    '1FC32410652A2C49BD845E567',    'FE0A9A7AFDAC44E42CB95D261', ...
    'B03062DC2B71995D5AD8B7DBE',    'F6C398993F598E2DF4235D3D5', ...
    '1BB2FB8B5BF24395C2EF3C5A1',    '2F920687D238CC7046EF6AFC9', ...
    '34163886FC4ED7F2A92EFDBB8',    '66A872CE47833FB2DFD5625AD', ...
    '99D5A70162C920A4BB9DE1CA8',    '81D71BD6E069A7ACCBEDC66CA', ...
    'A654524074A9E6780DB9D3EC6',    'C3396A101BEDAF623CFC5BB37', ...
    'C3D4AB211DF36F2111F2141CD',    '3DFF25EAE761739265AF145C1', ...
    '994909E0757D70CDE389102B5',    'B938535522D119F40C25FDAEC', ...
    'C71AB549C0491537026B390B7',    '0CDB8C9E7B53F55F5B0A0597B', ...
    '61C5FA252F1AF81144766494F',    '626027778FD3C6BB4BAA7A59D', ...
    'E745412FF53DEBD03F1C9A633',    '3592AC083F3175FA724639098', ...
    '52284D941C3DCAF2721DDB1FD',    '73B3D8F0AD55DF4FE814ED890', ...
    '94BF16C83BD7462F6498E0282',    'A8C3DE1AC668089B0B45B3579', ...
    'E23FFC2DD2C14388AD8D6BEC8',    'F2AC871CDF89DDC06B5960D2B', ...
    '06191EC1F622A77A526868BA1',    '22D6E2A768E5F35FFC8E01796', ...
    '25310A06675EB271F2A09EA1D',    '9F7993C621D4BEC81A0535703', ...
    'D62999EACF1C99083C0B4A417',    'F665A7EA441BAA4EA0D01078C', ...
    '46F3D3043F24CDEABD6F79543',    'E2E3E8254616BD96CEFCA651A', ...
    'E548231A82F9A01A19DB5E1B2',    '265C7F90A16F49EDE2AA706C8', ...
    '364A3A9EB0F0481DA0199D7EA',    '9810A7A898961263A0F749F56'];

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
E5aQ_secondary = [first_segement second_segement];

% Convert to bipolar format (0-->1; 1-->-1)
E5aQ_secondary = 1 - 2* E5aQ_secondary;

