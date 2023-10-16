function [F,L] = fusionframe(frametypes, varargin)
%FUSIONFRAME  Construct a new fusion frame
%   Usage: [F,L] = fusionframe(frametypes)
%          [F,L] = fusionframe(frametypes, 'weights',weights)
%          [F,L] = fusionframe(frametypes, 'ftype',params)
%
%   `F=fusionframe(ftype,...)` constructs a new fusion frame object *F* of type
%   *ftype*.
%
%   Input parameters:
%         frametypes     : Cell array with the desired frame types. For supported frames, see below.
%
%   Output parameters:
%         F              : fusion frame structure
%         L              : common frame length
%
%   `fusionframe` facilitates the construction of fusion frames by allowing
%   for their programmatic instantiation via the frame type of the local frames. 
%   Further specification of the frames is optionally ppossible. For a
%   list of supported frames as well as their defaults, see arg_fusionframe. 
%   Currently, only one set of optional parameters per fusion frame type is
%   possible.
%   The common frame length is calculated as the lowest common multiple of
%   the frame lengths of the local frames.
%
%   Additional parameters:
%
%     'weights',weights     the fusion frame weights; per default set to one
%
%     'frametype',ftype     further parameters for one or more specific
%                           local frames
%
%
%   Examples:
%   ---------
%
%   [F, L] = fusionframe({'dgt', 'fwt'}, 'dgt',{'gauss', 20, 30});
%
%
%   See also: frame


if ~iscell(frametypes)
   error('%s: The frame types need to be passed as a cell array.',upper(mfilename));
end

framenumber = numel(frametypes);

definput.keyvals.weights = ones(framenumber, 1);
definput.import={'fusionframe'};

[~,kv]=ltfatarghelper({},definput,varargin,'fusionframe');


supported_frames = fieldnames(kv);

fusionframe_args = [];
for ii = 1:framenumber
    for jj = 1:numel(supported_frames)
        frstr = [];
        if strcmp(frametypes{ii}, supported_frames{jj})
            input_args = kv.(supported_frames{jj});
            for kk = 1:numel(input_args)
                if ischar(input_args{kk})
                    frstr = strcat(frstr, sprintf(',''%s''',input_args{kk}));
                else
                    frstr = strcat(frstr, sprintf(',%i',input_args{kk}));
                end
            end
            %remove spurious first comma
            frstr = frstr(2:end);
            eval(sprintf('F_%i = frame(''%s'', %s);',ii, frametypes{ii},frstr));
            fusionframe_args = strcat(fusionframe_args, sprintf(", F_%i",ii));
            eval(sprintf('L(%i) = framelength(F_%i, 1);',ii, ii));
        end
    end
end

% calculate the common frame length
L_temp =zeros(framenumber, 1);
L_temp(1) = L(1);
for ll = 1:numel(L)-1
    L_temp(ll+1) = lcm(L_temp(ll), L(ll+1));
end
L = L_temp(end);

F = eval(sprintf("frame('fusion', kv.weights %s);",fusionframe_args));