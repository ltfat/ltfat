function [] = mulaclab(file,varargin)
%MULACLAB Graphical interface for audio processing using frame multipliers
%   Usage: mulaclab;
% 
%   When starting the interface, the user is asked to choose the processed
%   signal, named original signal in the interface.  Possible signals are
%   `.wav` files and `.mat` files containing decompositions preliminarily
%   saved using the `mulaclab` interface.  The interface only handles
%   monochannel signals. So for a multichannel `.wav` files, the first channel
%   is used as the original signal.
%
%   After choosing the original signal, the user is presented with the main
%   interface. This interface is divided in two areas: 
%
%   * The right part of the figure contains the visualizations, which
%     represent the spectrograms of the original and modified signals.
%
%     The 'Original signal' visualization is used to display the original
%     signal spectrogram and to graphically define the symbol of the Gabor 
%     multiplier that will be applied on this signal. 
%
%     The 'Overview of original signal' visualization also represents the
%     original signal spectrogram and can be used for fast zooming and moving
%     of the other visualizations. Zooming and moving is controlled by mouse
%     interaction with the white rectangle displayed on this visualization. 
%
%     The 'Modified signal' visualization is used to display the spectrogram
%     of the modified signal after application of the multiplier.
%
%     It is possible to hide the 'Overview of original signal' and 'Modified
%     signal' visulizations using the 'Visualization' menu.
%
%   * The left part of the figure contains panels with tools for user 
%     interaction.
%
%     The 'Audioplayer' panel contains the controls for audio playback of the
%     original and modified signal.
%
%     The 'Visualization' panel contains tools used to adapt the display
%     of the visualizations.
%
%     The 'Selection' panel contains tools and information concerning the
%     multilayered selection used to graphically specify the symbol of the 
%     multiplier.
% 
%   Known Matlab limitations: 
%
%   * When using Matlab on Linux with multiple screens, there might be a 
%     Matlab bug preventing the display of the multiplier symbol. This can be
%     solved by docking the figure.
%
%   * When using a Matlab version prior to 7.3 (R2006b), the rectangle
%     displayed on the 'Overview of original signal' visualization is not
%     automatically updated when using the zoom and pan tools of the 'Zoom' 
%     panel. It can be manually updated by re-clicking on the currently
%     selected tool or by changing the current tool.
%
%   The Matlab Image Processing Toolbox is required by the mulaclab function.
%
%   MULACLAB uses the GPC library available from
%   `<http://www.cs.man.ac.uk/~toby/gpc/>`_. This library is distributed
%   alongside LTFAT, but under different licensing conditions. Please see
%   the `ltfat/thirdparty/gpc/GPC-README.pdf` file for the exact conditions.

%
% questions : do we want to :
% - have the possibility to use Gaussian window? (any window?) - check
% - be able to export the symbol? 
% - be able to export the coeff? 
% - automatically generate a matlab file applying the corresponding gabor
% mutliplier? something like:
% [s, ...] = waveread('rightfile.wav');
% coeff = dgt(f, ...);
% load symbol.mat
% coeff = coeff.*symbol;
% res = idgt(coeff, ...);
%
% Known bugs:
% - visualization bug when the signal is zero everywhere
% - there might be a bug with noise level when using hole filling with
% noise, check what's happening. I tried to check, but I cannot find
% the bug. Still, something really strange is hapenning, for a given
% selection I will have a level which is too low, when everything works
% nicely with a slightly modified version of the selection.
% - when using undo it seems that parameters of the type of layer are not
% updated
%
% todo:
% - reorder the functions in this file to improve consistency
% - check the interface with matlab r14
% - check the looking of the interface under linux, windows and mac
% - add possibility to see individual atom of the frame (and the dual)
%   (with right clic menu?)
% - possibility to have a nice view (3D?) of the mask
% - possibility to choose the channel and a time selection of wave
% - switch between original and modified when playing
% - add editable keyboard shortcuts for almost everything
% - comment the functions
% - add a simple polygon simplification function (at least for straight 
% line with magic wand)
% - for changing mouse pointer, see what can be done using 
%   iptPointerManager
% - give possibility to use a Gaussian window
% - see if I change the scale used to display symbol (use dB?)
% - precise problems if public distribution of this file (use of 
%   PolygonClip and icons)
% - update information of contextual menu for layers to know what are the
%   currently selected properties
% - see if I put a limit in the number of levels of undo
% - add possibility to modify afterward polygons
%
% Acronym used in the comments:
% PI: possible improvement, used to specify things that could be done to
% improve the current version of mulaclab. In particular, it can be 
% efficiency improvements or suitable functionnality addition.
%
% _____________________ Description of shared data ________________________
%
% To share data between the different components of the Graphical User 
% Interface (GUI), nested functions are used
% (search for 'Sharing Data Among a GUI's Callbacks' in matlab's
% documentation for more information about this technique).
% So all the functionnalities of the GUI are implemented using nested
% functions, and to share data, a set of variables are defined at the top 
% level of the function (ie in the workspace of the function mulaclab).
%
% Here is a list of these variables with their description:
%
% sig: (struct) data describing the processed signal
%
% - sig.ori: (vector) samples of the original signal
%
% - sig.mod: (vector) samples of the modified signal
%
% - sig.sampFreq: (scalar) sampling frequency of the original and modified
%   signals
%
% - sig.nbBit: (scalar) number of bits to use when saving the modified
%   signal
%
% - sig.real: (bool) boolean to know if the signal is real (true)
%   or complex (false)
%
% frame: (struct) data describing the frame
%
% - frame.type: (string) type of the frame (currently only Gabor frames are 
%   implemented, other frames could be added, like wavelet frames,
%   nonstationnary Gabor frames, nonuniform filterbanks, general frames)
%
% - frame.def: (struct) data defining the frame (fields depend of the
%   frame type)
%
%   For Gabor frames, frame.def fields are:
%   - frame.def.winType: (string) type of window as defined in the firwin
%     function of LTFAT, see parameter 'name' in the help of firwin 
%     function for details
%
%   - frame.def.winLen: (scalar) window length in number of samples
%
%   - frame.def.hop: (scalar) hop size (ie Gabor frame time shift) in
%     number of samples (parameter 'a' of LTFAT dgt function)
%
%   - frame.def.nbFreq: (scalar) number of frequency bins (parameter 'M' 
%     of LTFAT dgt function)
%
% coeff: (struct) coefficients obtained by decomposition of the signal 
% on the frame
%
% - coeff.ori: (matrix, dimension 1: frequency, dimension 2:time) 
%   coefficients for the original signal sig.ori
%
% - coeff.mod: (matrix, dimension 1: frequency, dimension 2:time)
%   coefficients for the modified signal sig.mod
%
% - coeff.oriC: (matrix, dimension 1: frequency, dimension 2:time)
%   plotted spectrogram as image in DB
%
% - coeff.info: (struct)
%   additonal data necesary for the reconstruction
%
% sel: (struct) data describing the graphical selection
% The graphical selection contains several layers, each layer containing
% several polygons defining a region of the time-frequency that must be 
% modified. These polygons can be drawn freely by the user using the two 
% available tools (freehand selection and magic wand).
% sel contains data defining these polygons and layers, and data
% specifying how they should be plotted
%
% - sel.mode: (string) specifies how the next polygon drawn by the user 
%   must be combined with the polygons already defined in the currently
%   edited layer
%   possible modes are 'union', 'intersection', 'difference'
%
% - sel.curLay: index of the currently edited layer
%
% - sel.lay: (struct) data describing the layers, it's a vector of struct,
%   and sel.lay(ind) contains data for layer number ind
%
%   - sel.lay.convType: (string) conversion type of the layer, it specifies
%     how the values of the symbol must be computed for the current layer,
%     possible types are:
%
%     'constGain': apply a constant gain on the whole time-frequency region
%     corresponding to the polygons of the layer
%
%     'smoothBorder': also apply a gain on the region corresponding to the
%     polygons of the layer, but the gain values are smoothed at the border
%     to go linearly from 1 outside the polygons to the gain value
%     specified for the layer inside the polygons
%
%     'fill': try to fill the region corresponding to the polygons of the
%     layer according to the neighborhood level, this is done by
%     replacing the absolute value of the coefficient inside the polygons 
%     by the interpolated values using the level of the coefficient in a 
%     specified neigborhood outside the polygons
%     the phase of the coefficients are left unchanged
%
%     'fillNoise': same as 'fill', but the coefficients inside the polygons
%     of the layers are replaced by coefficients obtained from a noise and
%     multiplied by the interpolated absolute values
%
%   - sel.lay.param: (struct) parameters for the conversion of the layer,
%     it's a struct array, sel.lay(indLay).param(indParam) defines the
%     parameter number indParam of the layer number indLay
%     the number of parameters and their nature depends of the convType
%     parameter of the layer according to the following:
%
%     * if convType is 'constGain': only one parameter named 'Gain' is 
%       defined and it contain the value of the constant gain applied to 
%       region corresponding to the polygons of the layer
%
%     * if convType is 'smoothBorder': two parameters named 'Gain' and 
%       'Border' are defined
%       'Gain' contain the value of the gain applied inside region 
%       corresponding to the polygons of the layer away from the border
%       'Border' specify the distance from the border (in number of
%       coeffiicents) in which the gain value is linearly interpolated 
%       between 1 and the value given by 'Gain'
%
%     * if convType is 'fill' or 'fillNoise' : two parameters 'Height' and
%       'Width' are defined
%       these parameters specifiy the size of the ellipse used to define
%       the size of the neighborhood used the estimate the level around the
%       polygons
%       'Height' is the radius (in number of coefficients) of the ellipse
%       along frequency (y axis), and 'Width' is the radius (in number of 
%       coefficients) of the ellipse along time (x axis)
%       
%     - sel.lay.param.name: (string) name of the parameter, as displayed 
%       in the GUI under the layer list
%
%     - sel.lay.param.val: (scalar) the value of the parameter
%
%   - sel.lay.lineWidth: (scalar) width of the lines representing the 
%     polygons of the layer
%
%   - sel.lay.color: (colorspec) color of the lines representing the 
%     polygons of the layer
%
%   - sel.lay.marker: (string) marker type of the lines representing the 
%     polygons of the layer
%
%   - sel.lay.markerSize: (scalar) marker size of the lines representing 
%     the polygons of the layer
%
%   - sel.lay.lineStyle: (string) line style of the lines representing the 
%     polygons of the layer
%
%   - sel.lay.poly: (struct) data defining the polygons of the layer, it's
%     a struct array, sel.lay(indLay).poly(indPoly) defines the
%     polygon number indPoly of the layer number indLay
%
%     - sel.lay.poly.x: (row vector) x-coordinates (ie time-coordinates)
%       of the polygon points
%
%     - sel.lay.poly.y: (row vector) y-coordinates (ie 
%       frequency-coordinates) of the polygon points
%
%     - sel.lay.poly.hole: (bool) specifies if the polygon is a hole
%
%     - sel.lay.poly.id: (graphics object handle) handle to the line
%       representing the polygon
%
%   - sel.lay.label: (string) layer label as displayed in the layer list in
%     the GUI
%
%
% visu: (struct) data defining parameters of the visualizations
% visu is a struct array, visu(ind) contains data for the visualization
% number ind.
% The visualizations are the representations plotted on the right part of
% the GUI. They are plotted in the order specified by visu from top to
% bottom, that is to say that visu(1) is plotted on top, visu(2) is plotted
% under visu(1) and so on.
% PI: Currently only three visualizations are available: the spectrogram
% of the original signal (used as the support for graphical definition of 
% the multiplier), the spectrogram of the modified signal to observe the 
% effect of the multiplier, and an overview spectrogram of the original 
% signal for fast navigation in the time-frequency plane.
% But it should be easy to add some new visualizations if needed. To add a
% new visualization, adding an element to vector visu with appropriate 
% parameters and implementing the associated updating function should be
% enough. If the possibility to show and hide the new visualization is
% also needed, some handling of it in the menu is also required.
%
% - visu.label: (string) the title appearing on top of the visualization
%
% - visu.height: (scalar) number specifying the height of the
%   visualization. It is a number between 0 and 1 specifying the ratio of
%   the figure height that should be used for each visualisation. To cover
%   the whole figure, values should be chosen so that they sum up to 1 when
%   considering all the visible visualizations
%
% - visu.position: (vector) the position vector of the uipanels containing
%   the visualizations. It is automatically computed by function resize
%   using the values of visu.height.
%
% - visu.updateFcn: (function handle) the handle of the function that must
%   be used to update the plot of the visualization
%
% - visu.dynamic: (scalar) colorscale dynamic of the visualization
%   (positive number in dB)
%
% - visu.visible: (bool) boolean specifying if the visualization is
%   visible and must be plotted (true) or not (false)
%
% - visu.linkX: (bool) boolean specifying if the visualization x axe (time
%   axe) must be linked with the one of the other visualizations
%
% - visu.linkY: (bool) boolean specifying if the visualization y axe 
%   (frequency axe) must be linked with the one of the other visualizations
%
% - visu.linkCLim: (bool) boolean specifying if the visualization 
%   colorscale limits must be linked with the one of the other 
%   visualizations
%
% - visu.uipanelId: (graphics object handle) handle of the uipanel
%   containing the visualization
%
% - visu.axesId: (graphics object handle) handle of the axes of the
%   visualization
%
% gui: (struct) graphical user interface data 
% gui contains data used to define properties of certain graphics object
% (in particular their position) and data used to dynamically handle the
% inferface (in particular many graphics object handles)
%
% - gui.mainFigPos: (4x1 vector) position of the main figure using the
%   usual matlab position definition of the form [left bottom width height]
%
% - gui.fontSize: (scalar) font size (in points) for the elements appearing
%   in the panels on the left part of the interface
%
% - gui.textHeight: (scalar) height (in number of pixels) for text
%   uicontrol appearing in the panels on the left part of the interface
%
% - gui.textWidth: (scalar) width (in number of pixels) for text
%   uicontrol appearing in the panels on the left part of the interface
%
% - gui.editWidth: (scalar) width (in number of pixels) for edit
%   uicontrol appearing in the panels on the left part of the interface
%
% - gui.buttonBackgroundColor: (colorspec) background color used for
%   uicontrol appearing in the panels on the left part of the interface
%
% - gui.buttonWidth: (scalar) width (in number of pixels) for button
%   uicontrol appearing in the panels on the left part of the interface
%
% - gui.buttonHeight: (scalar) height (in number of pixels) for button
%   uicontrol appearing in the panels on the left part of the interface
%
% - gui.horiDist: (scalar) horizontal distance (in number of pixels) 
%   between two elements appearing in the panels on the left part of the 
%   interface
%
% - gui.vertDist: (scalar) vertical distance (in number of pixels) 
%   between two elements appearing in the panels on the left part of the 
%   interface
%
% - gui.margin: (1x2 vector) distance from the border of the panel (in 
%   number of pixels) of the element closest to the border of the panel.
%   First element is horizontal distance, second element is vertical
%   distance.
%
% - gui.marginSub: (1x2 vector) same as gui.margin but for elements 
%   appearing in sub-panels
%
% - gui.panelWidth: (scalar) width (in number of pixels) of the panels 
%   appearing on the left part of the interface
%
% - gui.panelHeight: (vector) height (in number of pixels) of the panels 
%   appearing on the left part of the interface. These values must be
%   adapted by hand if other parameters such as gui.buttonHeight are 
%   changed
%
% - gui.panelTitle: (cell array of strings) title of the panels appearing 
%   on the left part of the interface
%
% - gui.visuResizeButtonWidth: (scalar) width (in number of pixels) of 
%   the draggable resize buttons appearing between two visualizations on 
%   the right part of the interface
%
% - gui.visuResizeButtonHeight: (scalar) height (in number of pixels) of 
%   the draggable resize buttons appearing between two visualizations on 
%   the right part of the interface
%
% - gui.curTool: (string) currently used graphical tool, can be one of the
%   possible names stored in gui.tool.name, which are:
%   'free' for freehand selection
%   'level' for level selection (ie magic wand)
%   'zoomIn' for zoom in mode
%   'zoomOut' for zoom out mode
%   'pan' for pan mode
%
% - gui.tool: (struct) struct array containing tool data
%   gui.tool(ind) contains data for tool number ind
%
%   - gui.tool.buttonId: (graphics object handle) handle of the tool button
%
%   - gui.tool.name: (string) name of the tool, used to identify the
%     current tool in gui.curTool
%
%   - gui.tool.function (function handle) handle of the function that must
%     be executed to initialize the tool when clicking the corresponding 
%     button
%
%   - gui.tool.param: (struct) struct array containing the tool parameters
%     the content of this element depends on the tool and is empty for all
%     the tools except for the levele selection (magick wand), for which 
%     one parameter named 'Tolerance' is defined
%
%     - gui.tool.param.name: (string) the name of the parameter, as
%     displayed under the tools subpanel
%
%     - gui.tool.param.val: (string) the value of the parameter represented
%       as a string (so a conversion is needed if the parameter is a
%       number)
%
% - gui.symbOpacity: (scalar) opacity value used for plotting of symbol 
%   (any value between 0:transparent and 1:opaque)
%
% - gui.showSymb: (bool) boolean specifying if the symbol should be plotted
%   (true) or not (false)
%
% - gui.symbImageId: (graphics object handle) handle of the image used to
%   plot the symbol (the symbol is plotted as a transparent image on top 
%   of the image of the spectrogram)
%
% - gui.player: (struct) structure containing the handles associated with
%   the gui of the audioplayer
%
%   - gui.player.buttongroupId: (graphics object handle) handle of the
%     buttongroup used for the selection of the signal played (original or 
%     modified) 
%
%   - gui.player.buttonOriId: (graphics object handle) handle of the radio
%     button for the selection of original signal
%
%   - gui.player.buttonModId: (graphics object handle) handle of the radio
%     button for the selection of modified signal
%
% - gui.ResizeVisuButtonId: (graphics object handle vector) vector
%   containing the handles of the draggable resize buttons appearing 
%   between two visualizations on the right part of the interface
%
% - gui.mainFigId: (graphics object handle) handle of the main figure of 
%   the interface
%
% - gui.undoMenuId: (graphics object handle) handle of the 'Undo' menu
%   element
%
% - gui.redoMenuId: (graphics object handle) handle of the 'Redo' menu 
%   element
%
% - gui.visuMenu: (struct) structure containing menu handles used to show
%   and hide visualizations
%
%   - gui.visuMenu.showModId: (graphics object handle) handle of the 
%   'Show modified signal' menu element
%
%   - gui.visuMenu.showOverviewId: (graphics object handle) handle of the 
%   'Show overview of original signal' menu element
%
% - gui.exploreRect: (struct) handle of the graphical objects used to
%   represent the rectangular selection on the overview of original signal
%   visualization
%
%   - gui.exploreRect.patchId: (graphics object handle) handle of the patch
%   used to draw the movable selection rectangle
%
%   - gui.exploreRect.lineId: (graphics object handle) handle of the line
%   used to draw the movable points at the corners of the rectangle
%
% - gui.panelId: : (graphics object handle vector) vector containg the
%   handles of the uipanels drawn on the left part of the interface
%
% - gui.editDynamicId: (graphics object handle) handle of the edit
%   uicontrol used to specify the colorscale dynamic
%
% - gui.layListId: (graphics object handle) handle of the listbox
%   uicontrol used to show the slection layers list
%
% - gui.layPanelId: (graphics object handle) uipanel handle of the 'Layers'
%   subpanel
%
% - gui.toolPanelId: (graphics object handle) uipanel handle of the 'Tools'
%   subpanel
%    
% link: (struct) data for linking of axes of the different visualization
%
% - link.CLim: (linkprop object) link object used to link the colorscale 
%   limits (through the CLim property) of the visualizations having their
%   parameter visu.linkCLim set to true
%
% - link.XLim: (linkprop object) link object used to link the time axe 
%   (through the XLim property) of the visualizations having their 
%   parameter visu.linkXLim set to true
%
% - link.YLim: (linkprop object) link object used to link the frequency axe 
%   (through the YLim property) of the visualizations having their 
%   parameter visu.linkYLim set to true
%
% player: (struct) data used for audio playing
% 
% - player.ori: (audioplayer object) audioplayer used to play the original
%   signal
%
% - player.mod: (audioplayer object) audioplayer used to play the modified
%   signal
%
% - player.selected: (string) indicate the currently selected audio signal
%   can be 'ori' for the original signal, or 'mod' for the modified signal
%
% - player.loop: (bool) indicate if the audio signal must be looped (true)
%   or not (false) when playing
%
% - player.position: (scalar) sample position at which the audio playing
%   will start at next pressing of start
%
% - player.forceStop: (bool) indicate to the stop function of the 
%   audioplayer that the player must be stopped (needed to be able to stop
%   when looping)
%
% - player.positionPlotId: (graphics object handle vector) contains handles
%   to the moving vertical white lines that are displayed on the 
%   visualizations to show the current position when playing
%
% default default data for some variables
% PI: This could be used more systematically to define the default value
% of all the parameters of the interface. It could then be used to
% configure the application by loading variable default from a file during
% initialization, and we could also give the possibility to the user
% to modify these default values
%
% - default.dynamic: (scalar) default value for colorscale dynamic 
%   (positive number in dB)
%
% - default.opacity: (scalar) default value for opacity used for plotting
%   of symbol (any value between 0:transparent and 1:opaque)
%
% - default.sel: (struct) default sel variable (see the description of sel)
%
% - default.frame: (struct) default frame variable (see the description of
%   frame)     
%
% undoData: (cell) data for undo functionnality
% Each cell of undoData contain a copy of the sel variable as it was at a
% preceding state
%
% redoData: (cell) data for redo functionnality
% Each cell of redoData contain a copy of the sel variable as it was at a
% preceding state
%
% ___________________ End of description of shared data ___________________
%
% __________________________ Beginning of code ____________________________
%
% define the variables that will be used to share data between
% interface components. They are defined at the top level of the function
% so that they are usable by all the functions used by the interface (as
% these are nested functions).
coeff = struct;
default = struct;
frame = struct;
gui = struct;
link = struct;
player = struct;
redoData = {};
sel = struct;
sig = struct;
undoData = {};
visu = struct;
symbol = struct;
visucommon = struct;
export = struct;
export.xLim = 800;
export.limitXaxesRes = true;
export.symbol = {};

% Dictionary of structures field names string expansions
mulaclabDict = struct('winType','Window type',...
                      'nbFreq','Number of frequency bins',...
                      'hop','Hop size',...
                      'winLen','Window length',...
                      'wavelet','Wavelet type',...
                      'J','Decomposition depth');
                   
% Cellaray of supported frames                  
supportedFramesIdx = 1;
supportedFrames = ...
{...
  struct('type','Gabreal','def',...
    struct('winType','hann','winLen',256,'hop',round(256 / 4),...
           'nbFreq',256)),... 
  ...
  struct('type','Gabor','def',...
    struct('winType','hann','winLen',256,'hop',round(256 / 4),...
           'nbFreq',256)),...
  ...
  struct('type','DWT','def',...
    struct('wavelet','sym10','J',8)),...
  ...
  struct('type','UDWT','def',...
    struct('wavelet','sym10','J',8)),...
  ...
  struct('type','WFBT','def',...
    struct('wavelet','sym10','J',7)),... 
  ...
  struct('type','UWFBT','def',...
    struct('wavelet','sym10','J',4)),...   
};

iconpath=[ltfatbasepath,'mulaclab',filesep,'icons',filesep];

% initialize the interface
if(nargin<1)
   file = [];
end
initialize(file);

% end of main function, from here everything is handled with the callbacks
% using the following nested functions

% __________________________ INITIALIZATION _______________________________
  function initialize(file)

    % check that we have all the toolboxes and functions that we need
    checkToolbox;
    
    % Display spash screen
    warningText = ['\newline{\color{red}\bf PLEASE NOTE: When using a free-hand selection tool, do not \newline hold down the mouse button, but click on the corner points. }\newline'];
    msgboxText = [' MulacLab: Matlab GUI to manipulate the short-time Fourier ' ...
               'transform of a signal \newline using Gabor multipliers. '...
              '\newline MulacLab makes use of the GPC Polygon Clipping library ' ...
               'available from \newline http://www.cs.man.ac.uk/~toby/alan/software/'...
             'GPC and mulaclab are free for \newline non-commercial use only. Please see the ' ...
               'GPC homepage for the exact licensing \newline conditions.'];
    
    h=msgbox(msgboxText,'MulacLab licensing','modal');
    posh = get(h,'Position');
    posh(end) = 120;
    set(h,'Position',posh);
    set(findall(h,'type','text'), 'Interpreter','tex', 'string', [msgboxText,warningText]); 
    uiwait(h,16);
    if ishandle(h)
       % Timeout occured, we close the box automatically.
      close(h);
    end;
    
    loadFile = 0;
    if nargin==0 || isempty(file)
         % get an original signal or decomposition to have something to show
          [fileName, pathName, filterIndex] = uigetfile('*.wav;*.mat',...
            'Open original signal or decomposition'); 

          if fileName == 0
            % the user pressed cancel, we stop here
            return;
          end
          loadFile = 1;
     elseif(ischar(file))
          fileName = file;
          pathName = '';
          loadFile = 1;
    end     

    if loadFile
      [~, ~, ext] = fileparts(fileName);
      switch lower(ext)
        case '.wav'
          % read the orginal signal in a wave file
          [sig.ori, sig.sampFreq, sig.nbBit] =...
            wavread([pathName, fileName]);
          sig.real = true;
        case '.mat'
           % read the original signal decomposition in a mat file
          data = load([pathName, fileName]);
          if isfield(data, 'frame') && isfield(data, 'savedSig')
           frame = data.frame;
           sig = data.savedSig;
          else
           errordlg(...
             [fileName ' is not a valid signal decomposition file']);
           return;
          end
      end
    elseif(nargin>=1)
       definput.keyvals.Fs=[];
       [flags,kv,Fs]=ltfatarghelper({'Fs'},definput,varargin);
       if(isempty(Fs))
          error('%s: Second parameter: sampling frequency is missing. ',upper(mfilename));
       end
       
       if(numel(find(size(file)>1))>1)
          error('%s: Input has to be one channel signal.',upper(mfilename));
       end
       %normalize and change to column vector
       sig.ori=file(:)/max(abs(file));
       sig.real=isreal(file);
       sig.sampFreq = Fs;
       sig.nbBit = 16;
    else
        errordlg(sprintf('%s: Unrecognized input parameters.',upper(mfilename))   );
           return;
    end

    if size(sig.ori, 2) > 1
       % multichannel wave, we keep only the first channel
       sig.ori = sig.ori(:, 1);
    end
    
    sig.mod = sig.ori;
    % intialize the different parameters of the application
    initializeParameter;
    
    % initialize the frame definition and the coefficients
    frame = default.frame;
    %coeff.ori = calculateCoeff(sig.ori);
    %coeff.mod = coeff.ori;
    
    % intialize audioplayers
    %initializePlayer;
    
    % create main window
    gui.mainFigId = figure(...
      'Name', 'MulAcLab',...
      'Position', gui.mainFigPos,...
      'Toolbar', 'none',...
      'Colormap', jet(256),...
      'ResizeFcn', @resize,...
      'MenuBar', 'none',...
      'WindowStyle', 'normal',...
      'Visible', 'off');

    % create menu of main window
    menuFileId = uimenu(gui.mainFigId, 'Label','File');
    uimenu(menuFileId,...
      'Label','Open original signal',...
      'Callback',@openOriSig);
    uimenu(menuFileId,...
      'Label','Import original signal decomposition',...
      'Callback',@openOriDecompo);
    uimenu(menuFileId,...
      'Label','Import selection',...
      'Callback',@importSel);
    uimenu(menuFileId,...
      'Label','Import symbol',...
      'Callback',@importSymbol);
    uimenu(menuFileId,...
      'Label','Save modified signal',...
      'Callback',@saveModSig,...
      'Separator', 'on');
    uimenu(menuFileId,...
      'Label','Save original signal decomposition',...
      'Callback',@saveOriDecompo);
    uimenu(menuFileId,...
      'Label','Save modified signal decomposition',...
      'Callback',@saveModDecompo);
    uimenu(menuFileId,...
      'Label','Export selection as ...',...
      'Callback',@exportSel);
    uimenu(menuFileId,...
      'Label','Export symbol as ...',...
      'Callback',@exportSymbol);

    menuEditId = uimenu(gui.mainFigId,...
      'Label','Edit');
    gui.undoMenuId = uimenu(menuEditId,...
      'Label','Undo',...
      'Accelerator', 'z',...
      'enable', 'off',...
      'Callback',@undo);
    gui.redoMenuId = uimenu(menuEditId,...
      'Label','Redo',...
      'Accelerator', 'y',...
      'enable', 'off',...
      'Callback',@redo);

    menuFrameId = uimenu(gui.mainFigId,...
      'Label','Frame');
    menuFrameTypeId = uimenu(menuFrameId,...
      'Label','Choose frame type');
    for ind=1:length(supportedFrames)
       uimenu(menuFrameTypeId,...
      'Label',supportedFrames{ind}.type,...
      'Callback',@changeFrameType,...
      'UserData',ind);
    end
    uimenu(menuFrameId,...
      'Label','Edit frame parameters',...
      'Callback',@changeFrameDef);
    uimenu(menuFrameId,...
      'Label','Load frame parameters',...
      'Callback',@loadFrameParam,...
      'Separator', 'on');
    uimenu(menuFrameId,...
      'Label','Save frame parameters',...
      'Callback',@saveFrameParam);

    menuVisuId = uimenu(gui.mainFigId,...
      'Label','Visualization');
    % !!! the folowing checked on or off should be set automatically
    % during initialization
    gui.visuMenu.showModId = uimenu(menuVisuId,...
      'Label','Show modified signal',...
      'Checked', 'off',...
      'Callback',@toggleModVisu);
    gui.visuMenu.showOverviewId = uimenu(menuVisuId,...
      'Label','Show overview of original signal',...
      'Checked', 'on',...
      'Callback',@toggleOverviewVisu);
     
    % intialize visualization
    %updateVisu(true);

    % create gui panels
    currentY = 1;
    for indPanel = 1:length(gui.panelHeight)
      gui.panelId(indPanel) = uipanel(gui.mainFigId, ...
        'Title', gui.panelTitle{indPanel},...
        'TitlePosition', 'centertop',...
        'FontSize', gui.fontSize, ...
        'Units', 'pixels', ...
        'Position', [1 currentY gui.panelWidth gui.panelHeight(indPanel)]);
      currentY = currentY + gui.panelHeight(indPanel);
    end

    % create audioplayers panel
     drawAudioplayer(gui.panelId(1));
% 
% 
%     % create visualization tools panel
     drawVisualizationTool(gui.panelId(2));
% 
%     % create selection tools panel
     drawSelectionTool(gui.panelId(3));

    % activate the current tool
    % changeTool([], [], gui.curTool);

    % end of initialization, show the interface to the user
    set(gui.mainFigId, 'Visible', 'on');
    resetMulaclab();
    resetSel();
    resetSymbol();
  end

  function initializeParameter()
    % set default values

    % default colorscale dynamic in dB
    default.dynamic = 60;
    
    % default opacity parameter for symbol plotting (from 0:transparent to
    % 1:opaque)
    default.opacity = 0.5;
    
    % default selection
    default.sel.mode = 'union'; 
    default.sel.curLay = 1;
    default.sel.lay.convType = 'constGain';
    default.sel.lay.param(1).name = 'Gain';
    default.sel.lay.param(1).val = 0;
    default.sel.lay.lineWidth = 2;
    default.sel.lay.color = 'w';
    default.sel.lay.marker = 'none';
    default.sel.lay.markerSize = 2;
    default.sel.lay.lineStyle = '-';
    default.sel.lay.poly = [];
    default.sel.lay.label = 'Layer';
    default.frame = supportedFrames{supportedFramesIdx};
    
    default.symbol.data.name = 'Selection';
    default.symbol.data.val = [];
    default.symbol.data.invert = false;
    default.symbol.curSymb = 1;

    
    default.zerotodb = -40;
    
    % define parameters for the gui (graphical user interface)

    monitorPos = get(0,'monitorposition');
    monitorPos = monitorPos(1, :);
    ratio = 0.75;
    gui.mainFigPos = round([monitorPos([3, 4])*0.5*(1-ratio),...
      round(monitorPos([3, 4])*ratio)]);
    gui.fontSize = 9;
    gui.textHeight = gui.fontSize + 7;
    gui.textWidth = 57;
    gui.editWidth = 28;
    gui.buttonBackgroundColor = [0.7 0.7 0.7]; % chosen to fit icons color
    gui.buttonWidth = 27;
    gui.buttonHeight = 22;
    gui.horiDist = 3;
    gui.vertDist = 2;
    gui.margin = [5, 5];
    gui.marginSub= [3, 3];
    nbButtton = 3; % number of buttons on one line
    gui.panelWidth = 2*gui.margin(2) + (nbButtton-1) * gui.horiDist +...
      nbButtton*gui.buttonWidth;

    gui.panelHeight(1) = 91;
    gui.panelTitle{1} = 'Audioplayer';

    gui.panelHeight(2) = 220;
    gui.panelTitle{2} = 'Visualization';

    gui.panelHeight(3) = 260;
    gui.panelTitle{3} = 'Selection';

    gui.visuResizeButtonWidth = 15;
    gui.visuResizeButtonHeight = 15;

    gui.curTool = 'freehand';

    gui.tool = struct;
    gui.symbOpacity = 0.5;
    gui.showSymb = false;
    gui.symbImageId = [];

    gui.player.buttongroupId = [];
    
    % define parameters for the visualizations

    % !!! height of visu at initilisation must be choosen so that the sum
    % of visible visu is 1
    
    % Main visualization: spectrogram of the original signal
    visu(1).label = 'originalMain';
    visu(1).height = 0.67;
    visu(1).position = [0, 0, 1, 1]; % we can initialise this with wathever
                                     % value as it will be corrrected 
                                     % during resize
    visu(1).updateFcn = @updateVisuOriSpec;
    visu(1).dynamic = default.dynamic;
    visu(1).visible = true; % Important: the orginalMain visu must always
                            % stay visible
    visu(1).linkX = true;
    visu(1).linkY = true;
    visu(1).linkCLim = true;
    visu(1).uipanelId = [];
    visu(1).axesId = [];

    % Visualization of the spectrogram of the modified signal
    visu(2).label = 'modified';
    visu(2).height = 0.;
    visu(2).position = [0, 0, 1, 1];
    visu(2).updateFcn = @updateVisuModSpec;
    visu(2).dynamic = default.dynamic; % useless if linkCLim is true
    visu(2).visible = false;
    visu(2).linkX = true;
    visu(2).linkY = true;
    visu(2).linkCLim = true;
    visu(2).uipanelId = [];
    visu(2).axesId = [];

    % Visualization for the overview of the original spectrogram with fast
    % scrolling
    visu(3).label = 'originalOverview';
    visu(3).height = 0.33;
    visu(3).position = [0, 0, 1, 1];
    visu(3).updateFcn = @updateVisuOriOverview;
    visu(3).dynamic = default.dynamic; % useless if linkCLim is true
    visu(3).visible = true;
    visu(3).linkX = false;
    visu(3).linkY = false;
    visu(3).linkCLim = true;
    visu(3).uipanelId = [];
    visu(3).axesId = [];
    
    gui.ResizeVisuButtonId = NaN(length(visu), 1);
  end

  function checkToolbox()

    if isoctave
        error('This function requires Matlab.');      
    end;
    
    % Check that the Image Processing Toolbox is available
    if isempty(ver('images'))
        error(...
          'This function requires the Matlab Image Processing Toolbox.');
    end

    % Check that the function for polygon clipping is available
    temp1.x = [1];
    temp1.y = [1];
    temp1.hole = false;
    temp2 = temp1;
    try
      temp = PolygonClip(temp1, temp2);
    catch 
      error(['This function requires the function PolygonClip for ' ...
             'polygon clipping. Please compile using the "ltfatmex" ' ...
             'command.']);
    end

    % !!! check if some other toolboxes are needed (signal processing?)
  end

% ___________________ COMPUTATIONAL FUNCTIONS _____________________________

  function coef = calculateCoeff(f)
     switch lower(frame.type)
       case 'gabreal'
          win = firwin(frame.def.winType,frame.def.winLen);
          coef = dgtreal(f, win, frame.def.hop, frame.def.nbFreq);
       case 'gabor'
          win = firwin(frame.def.winType,frame.def.winLen);
          coef = dgt(f, win, frame.def.hop, frame.def.nbFreq);
       case 'dwt'
          [coef, coeff.info] = fwt(f,frame.def.wavelet,frame.def.J);
       case 'wfbt'
          [coef, coeff.info] = wfbt(f,{frame.def.wavelet,frame.def.J,'full'});
       case 'uwfbt'
          [coef, coeff.info] = uwfbt(f,{frame.def.wavelet,frame.def.J,'full'});
       case 'udwt'
          [coef, coeff.info] = ufwt(f,frame.def.wavelet,frame.def.J);
       otherwise
          error('%s: Unrecognized frame type',upper(mfilename));
     end
  end

  function fhat = recFromCoeff(coef)
     switch lower(frame.type)
       case 'gabreal'
          win = gabdual(firwin(frame.def.winType,frame.def.winLen), ...
          frame.def.hop, frame.def.nbFreq);
          fhat = idgtreal(coef, win, frame.def.hop,frame.def.nbFreq, length(sig.ori));
       case 'gabor'
          win = gabdual(firwin(frame.def.winType,frame.def.winLen), ...
          frame.def.hop, frame.def.nbFreq);
          fhat = idgt(coef, win, frame.def.hop,length(sig.ori));
       case 'dwt'
          fhat = ifwt(coef,coeff.info);
       case 'wfbt'
          fhat = iwfbt(coef,{frame.def.wavelet,frame.def.J,'full'},length(sig.ori));
       case 'uwfbt'
          fhat = iuwfbt(coef,{frame.def.wavelet,frame.def.J,'full'});
       case 'udwt'
          fhat = iufwt(coef,coeff.info);
       otherwise
          error('%s: Unrecognized frame type',upper(mfilename)); 
     end
    fhat = real(fhat);
  end

  function C = plotCoeff(coef,axesId,key,value)
    switch lower(frame.type)
       case 'gabreal'
          C = plotdgtreal(coef,frame.def.hop,frame.def.nbFreq,sig.sampFreq,key,value);
       case 'gabor'
          C = plotdgt(coef,frame.def.hop,sig.sampFreq,key,value);
       case {'dwt','udwt','wfbt','uwfbt','ufwt'}
          C = plotwavelets(coef,coeff.info,sig.sampFreq,key,value);
       otherwise
          error('%s: Unrecognized frame type',upper(mfilename));
    end
    if export.limitXaxesRes
        C=interp1(linspace(0,1,size(C,2)),C.',linspace(0,1,export.xLim),'nearest');
        C=C.';
    end
  end

  function mult = convSymbToCoefFormat(symb)
     switch lower(frame.type)
        % Classical gabor coeff format 
        case {'gabreal','gabor'}
            if export.limitXaxesRes
              L = size(coeff.ori,2);
              Lplot = size(symb,2);
              mult=interp1(linspace(0,1,Lplot),symb.',linspace(0,1,L),'nearest');
              mult = mult.';
            else
              mult = symb;
           end
           
        % Cell array containing column vectors format
        case 'wfbt'
           M = numel(coeff.ori);
           Lplot = size(symb,2);
           Lc = cellfun(@(cEl) size(cEl,1),coeff.ori);
           mult = cell(size(coeff.ori));
           for m=1:M
              mult{m} = interp1(linspace(0,1,Lplot),symb(m,:),linspace(0,1,Lc(m)),'nearest').';
           end
        % FWT specific coefficient format
        case 'dwt'
           Lc = coeff.info.Lc;
           Lcstart = cumsum([1;Lc(1:end-1)]);
           Lcend = cumsum(Lc);
           mult = zeros(sum(Lc),1);
           M = numel(Lc);
           Lplot = size(symb,2);
           for m=1:M
              mult(Lcstart(m):Lcend(m)) = interp1(1:Lplot,symb(m,:),linspace(1,Lplot,Lc(m)),'nearest').';
           end  
        % ufilterbak and all wavelet functions beginning with u
        case {'udwt','uwfbt'}
           M = size(coeff.ori,2);
           L = size(coeff.ori,1);
           Lplot = size(symb,2);
           mult = zeros(size(coeff.ori));
           for m=1:M
              mult(:,m) = interp1(1:Lplot,symb(m,:),linspace(1,Lplot,L),'nearest').';
           end
        otherwise
          error('%s: Unrecognized frame type',upper(mfilename));
     end

  end

  function symb = convCoefFormatToSymb(mult)
     switch lower(frame.type)
        % Classical gabor coeff format 
        case {'gabreal','gabor'}
            if export.limitXaxesRes
              L = size(coeff.ori,2);
              Lmult = size(mult,2);
              symb=interp1(linspace(0,1,Lmult),mult.',linspace(0,1,L),'nearest');
              symb = symb.';
            else
              symb=mult;
           end

        % Cell array containing column vectors format
        case 'wfbt'
           M = numel(coeff.ori);
           Lsymb = size(coeff.oriC,2);
           Lc = cellfun(@(cEl) size(cEl,1),mult);
           symb = zeros(size(coeff.oriC));
           for m=1:M
              symb(m,:) = interp1(linspace(0,1,Lc(m)),mult{m},linspace(0,1,Lsymb),'nearest');
           end
        % FWT specific coefficient format
        case 'dwt'
           Lc = coeff.info.Lc;
           Lcstart = cumsum([1;Lc(1:end-1)]);
           Lcend = cumsum(Lc);
           symb = zeros(size(coeff.oriC));
           %mult = zeros(sum(Lc),1);
           M = numel(Lc);
           Lplot = size(symb,2);
           for m=1:M
             symb(m,:) = interp1(linspace(1,Lplot,Lc(m)),mult(Lcstart(m):Lcend(m)),1:Lplot,'nearest');
           end  
        % ufilterbak and all wavelet functions beginning with u
        case {'udwt','uwfbt'}
           L = size(coeff.oriC,2);
           M = size(coeff.ori,2);
           Lmult = size(mult,1);
           symb = zeros(size(coeff.oriC));
           for m=1:M
              symb(m,:) = interp1(1:Lmult,mult(:,m),linspace(1,Lmult,L),'nearest');
           end
        otherwise
          error('%s: Unrecognized frame type',upper(mfilename));
     end

  end

  function applyMultiplier(c,mult)
     if(isnumeric(c)&&isnumeric(mult))
        coeff.mod = c.*mult;
     elseif(iscell(c)&&iscell(mult))
        coeff.mod = cellfun(@(x,y) x.*y,c,mult,'UniformOutput',false);
     else
        error('%s: Unrecognized frame type',upper(mfilename));
     end
     
     sig.mod  = recFromCoeff(coeff.mod);    
     updatePlayerMod;
     coeff.mod = calculateCoeff(sig.mod);
     modInd = find(strcmp('modified', {visu.label}), 1);
     if ~visu(modInd).visible
       toggleModVisu;
     end
     updateVisu;
  end

  function applySel(objId, eventData)    
    % convert the selection to get the symbol of the multiplier and apply 
    % the multiplier
    if isempty(symbol.data(symbol.curSymb).val)    
      mult = convSymbToCoefFormat(convSelToSymb());
    else
      mult = convSymbToCoefFormat(symbol.data(symbol.curSymb).val);
    end
     applyMultiplier(coeff.ori, mult);
  end

  function symb = convSelToSymb()
     symb = ones(size(coeff.oriC));
    
    for indLay = 1:length(sel.lay)
      if ~isempty(sel.lay(indLay).poly)
        mask = convPolyToMask(sel.lay(indLay).poly);
        symbCurr = convMaskToSymb(mask, indLay);
        symb = symb.*symbCurr;
      else
        warning(['Selection layer ' num2str(indLay)...
          ' is currently empty']);
      end
    end

  end

  function mask = convPolyToMask(poly)
    mask = zeros(size(coeff.oriC)); 
    % !!! Note: type of data for the mask might be optmized (bool, uint8?)

    % !!! could be optimized by doing the conversion in smaller boxes 
    % around each polygon (not on the whole time-frequency plane)
    % and also working with a smaller mask when signal is real
    for ind = 1:length(poly)
      if poly(ind).hole
        mask = mask - poly2mask(round(convAxesToIndX(poly(ind).x)),...
          round(convAxesToIndY(poly(ind).y)),...
          size(mask,1), size(mask,2));
      else
        mask = mask + poly2mask(round(convAxesToIndX(poly(ind).x)),...
          round(convAxesToIndY(poly(ind).y)),...
          size(mask,1), size(mask,2));
      end
    end

    % !!! modifiy this: don't do the symetrisation on the mask, but on the
    % symbol to insure that the results is real
%     if sig.real
%       temp = ceil(size(mask, 1)/2);
%       mask(end:-1:end-temp+2,:) = mask(2:temp,:);
%     end
  end

  function symb = convMaskToSymb(mask, indLay)
    % !!! precise that a case should be added here when creating new
    % selection type
    switch sel.lay(indLay).convType
      case 'constGain'
        gain = sel.lay(indLay).param(1).val;
        symb = ones(size(mask));
        symb = symb + (gain-1)*mask;
      case 'smoothBorder'
        gain = sel.lay(indLay).param(1).val;
        threshold = sel.lay(indLay).param(2).val;
        % bwdist returns single
        symb = double(bwdist(~mask));
        symb(symb > threshold) = threshold;
        symb = symb * (gain-1)/threshold + 1;
      case 'fill'
        % !!! first filling test, just interpolating in 3D the modulus
        
    
        % estimation of the level in the neighbourhood of the selection
        
        timeRadius = sel.lay(indLay).param(1).val;
        freqRadius = sel.lay(indLay).param(2).val;
        
        % find the points that are just next to the selection
        neighbour = imdilate(mask, strel('square',3)) - mask;
        
        % for these points, compute mean level in a specified neighbourhood
        % (but not taking into account the values of the transform inside 
        % the selection)
        
        % create the matrix defining the neighbourhood
        mat = computeEllipse(freqRadius, timeRadius); 
        % !!! this could be a more general filter if we want
        %mat  = mat./norm(mat);
        %temp1 = roifilt2(mat, abs(coeff.ori).^2 .* (1-mask), neighbour);
        temp1 = roifilt2(mat, 10.^(coeff.oriC/20) .* (1-mask), neighbour);
        temp2 = roifilt2(mat, 1-mask, neighbour);
        
        % !!! specifiy in comment that we have to do this to take into
        % account the hole
        
        [ind1Neigh, ind2Neigh] = find(neighbour);

        valNeigh = sqrt(temp1(sub2ind(size(temp1), ind1Neigh,...
          ind2Neigh))./temp2(sub2ind(size(temp2), ind1Neigh, ind2Neigh)));

        [ind1Mask, ind2Mask] = find(mask);

        try 
          valInterp = griddata(ind1Neigh, ind2Neigh, valNeigh,...
            ind1Mask, ind2Mask, 'linear');
        catch
          % sometime griddata fail, and adding this option seem to solve
          % the problem, but I don't really know what I'm doing here !!!
          valInterp = griddata(ind1Neigh, ind2Neigh, valNeigh,...
            ind1Mask, ind2Mask, 'linear', {'QJ'});
        end

        symb = ones(size(mask));
        % !!! there will be division by zero problem in the following
        % solve this problem
        symb(sub2ind(size(symb), ind1Mask, ind2Mask)) = valInterp ./...
          abs(coeff.oriC(sub2ind(size(coeff.oriC), ind1Mask, ind2Mask)));

        % !!! temporary modif: should be removed when modification are done
        % so that symetrisation is done on symbol and not mask when signal
        % is real
%         if sig.real
%           temp = ceil(size(symb, 1)/2);
%           symb(end:-1:end-temp+2,:) = conj(symb(2:temp,:));
%         end
      case 'fillNoise'
        % other hole filling test, using coefficients taken from a noise
        % estimation of the level in the neighbourhood of the selection
        
        timeRadius = sel.lay(indLay).param(1).val;
        freqRadius = sel.lay(indLay).param(2).val;
        
        % find the points that are just next to the selection
        neighbour = imdilate(mask, strel('square',3)) - mask;
        
        % for these points, compute mean level in a specified neighbourhood
        % (but not taking into account the values of the transform inside 
        % the selection)
        
        % create the matrix defining the neighbourhood
        mat = computeEllipse(freqRadius, timeRadius);
        % !!! this could be a more general filter if we want
        
        temp1 = roifilt2(mat, 10.^(coeff.oriC/20).* (1-mask), neighbour);
        temp2 = roifilt2(mat, 1-mask, neighbour);
        
        % !!! specifiy in comment that we have to do this to take into
        % account the hole
        
        [ind1Neigh, ind2Neigh] = find(neighbour);

        valNeigh = sqrt(temp1(sub2ind(size(temp1), ind1Neigh,...
          ind2Neigh))./temp2(sub2ind(size(temp2), ind1Neigh, ind2Neigh)));

        [ind1Mask, ind2Mask] = find(mask);

        try 
          valInterp = griddata(ind1Neigh, ind2Neigh, valNeigh, ...
            ind1Mask, ind2Mask, 'linear');
        catch
          % sometime griddata fail, and adding this option seem to solve
          % the problem, but I don't really know what I'm doing here !!!
          valInterp = griddata(ind1Neigh, ind2Neigh, valNeigh,...
            ind1Mask, ind2Mask, 'linear', {'QJ'});
          
          % !!! remove this
          disp('Problem with griddata')
        end

        symb = ones(size(mask));

        % compute the transform of a noise
        % !!! should be done in a nicer way with an external function 
        % so that it works with any kind of frame
        %noiseLen = (max(ind2Mask) - min(ind2Mask) + 1) * ( numel(sig.ori)/size(coeff.oriC,2)) ;
        noise = rand(numel(sig.ori), 1) - 0.5;

        % win = firwin(frame.def.winType, frame.def.winLen); 
        % coeffNoise = frame.funcAna(noise, win, frame.def.hop, frame.def.nbFreq);
        coeffNoise = calculateCoeff(noise);
        % !!! test figure(100); imagesc(20*log10(abs(coeffNoise)))
        
        % !!! there will be division by zero problem in the following
        % solve this problem
        
        % !!! old version
%         symb(sub2ind(size(symb), ind1Mask, ind2Mask)) = valInterp .*...
%           coeffNoise(sub2ind(size(coeffNoise), ind1Mask,...
%           ind2Mask-min(ind2Mask)+1)) ./...
%           coeff.ori(sub2ind(size(coeff.ori), ind1Mask, ind2Mask)) ./...
%           abs(coeffNoise(sub2ind(size(coeffNoise), ind1Mask,...
%           ind2Mask-min(ind2Mask)+1)));



        % !!! here all the interpolations and mean are done directly on the
        % absolute value of the  coefficients, it could be more meaningfull
        % to work with the sqaure of the modulus due to its possible
        % intepretation a energy ditribution
        coeffNoiseSymb = convCoefFormatToSymb(coeffNoise);
        meanLevelNoise = sum(abs(coeffNoiseSymb(sub2ind(size(coeffNoiseSymb), ...
          ind1Mask, ind2Mask-min(ind2Mask)+1)))) / length(ind1Mask);
        
        % !!! tests
        
        
%         symb(sub2ind(size(symb), ind1Mask, ind2Mask)) = valInterp .*...
%           coeffNoise(sub2ind(size(coeffNoise), ind1Mask,...
%           ind2Mask-min(ind2Mask)+1)) ./...
%           (coeff.ori(sub2ind(size(coeff.ori), ind1Mask, ind2Mask)) *...
%           meanLevelNoise);

        symb(sub2ind(size(symb), ind1Mask, ind2Mask)) = valInterp .*...
          coeffNoiseSymb(sub2ind(size(coeffNoiseSymb), ind1Mask,...
          ind2Mask-min(ind2Mask)+1)) ./...
          (coeffNoiseSymb(sub2ind(size(coeffNoiseSymb), ind1Mask, ind2Mask)) *...
          meanLevelNoise);

%         meanLevelOri = sum(abs(coeff.ori(sub2ind(size(coeff.ori), ...
%           ind1Mask, ind2Mask)))) / length(ind1Mask);
% 
%         symb(sub2ind(size(symb), ind1Mask, ind2Mask)) = ...
%           meanLevelOri*coeffNoise(sub2ind(size(coeffNoise), ind1Mask,...
%           ind2Mask-min(ind2Mask)+1+ 5)) ./...
%           (coeff.ori(sub2ind(size(coeff.ori), ind1Mask, ind2Mask))*meanLevelNoise);

        % !!! temporary modif: should be removed when modification are done
        % so that symetrisation is done on symbol and not mask when signal 
        % is real
%         if sig.real
%           temp = ceil(size(symb, 1)/2);
%           symb(end:-1:end-temp+2,:) = conj(symb(2:temp,:));
%         end
        
    end % of switch
    
  end

  function xData = convIndToAxesX(ind)
    % convert scale from image index to scale used by axes along X axe
    % xData = (ind-1) * (frame.def.hop / sig.sampFreq);
    hop = numel(sig.ori)/size(coeff.oriC,2);
    xData = (ind) * (hop / sig.sampFreq);
  end

  function ind = convAxesToIndX(xData)
    % !!! there is no round (on purpose), should be added if we really want
    % an integer index
    % convert scale from scale used by axes to image index along X axe
    % ind = xData * (sig.sampFreq / frame.def.hop) + 1;
    hop = numel(sig.ori)/size(coeff.oriC,2);
    ind = xData * (sig.sampFreq / hop);
  end

  function yData = convIndToAxesY(ind)
    % convert scale from image index to scale used by axes along Y axe
    % yData = (ind-1) * (sig.sampFreq / frame.def.nbFreq);
    Ylim = visucommon.oriYlim;
    yData = (ind) * ((Ylim(2)-Ylim(1)) / (size(coeff.oriC,1))) +Ylim(1);
  end

  function ind = convAxesToIndY(yData)
    % !!! there is no round (on purpose), should be added if we really want
    % an integer index
    % convert scale from scale used by axes to image index along Y axe
    Ylim = visucommon.oriYlim;
    ind = (yData-Ylim(1)) * (size(coeff.oriC,1) / (Ylim(2)-Ylim(1)) );
  end

  function ellipse =computeEllipse(dim1Radius, dim2Radius)
    ellipse = zeros(2*dim1Radius+1, 2*dim2Radius+1);
    ellipse(dim1Radius+1, :) = 1;
    ellipse(:, dim2Radius+1) = 1;
    for ind1 = 1:dim1Radius-1
      for ind2 = 1:dim2Radius-1
        if ((ind1/dim1Radius)^2 + (ind2/dim2Radius)^2) <= 1
          ellipse(dim1Radius+1+ind1, dim2Radius+1+ind2) = 1;
          ellipse(dim1Radius+1-ind1, dim2Radius+1+ind2) = 1;
          ellipse(dim1Radius+1+ind1, dim2Radius+1-ind2) = 1;
          ellipse(dim1Radius+1-ind1, dim2Radius+1-ind2) = 1;
        end
      end
    end
  end
% _____________________ GENERAL FUNCTIONS FOR UI __________________________

  function matlabVersion = getMatlabVersion()
    matlabVersion = ver('Matlab');
    matlabVersion = matlabVersion.Version;
    % the version number is sometime of the form X.X.X, so it cannot be
    % easily converted into a number, so we need to modify it so that it 
    % is more usable for comparison
    temp = find(matlabVersion == '.');
    temp = temp(2:end);
    for  ind = 1:size(temp, 2)
      % we remove the extra points '.' in the string
      matlabVersion = matlabVersion([1:temp(ind)-1, temp(ind)+1:end]);
    end
    matlabVersion = str2num(matlabVersion);
  end

  function resetMulaclab()
    initializePlayer;
    coeff.ori = calculateCoeff(sig.ori);
    sig.mod = sig.ori;
    coeff.mod = coeff.ori;
    resetSel();
    undoData = {};
    modInd = find(strcmp('modified', {visu.label}), 1);
    if visu(modInd).visible
       toggleModVisu();
    end
    updateVisu(true); % update visualization with reset of axes limits
    
    oriInd = find(strcmp('originalMain', {visu.label}), 1);
    visucommon.oriYlim = get(visu(oriInd).axesId,'Ylim');
    
    % activate the current tool
    changeTool([], [], gui.curTool);
    
  end

  function openOriSig(objId, eventData)   
    [fileName, pathName, filterIndex] = uigetfile('*.wav',...
      'Open original signal');
    if fileName == 0
      % user pressed cancel, stop here
      return;
    end
    [newSig, sampFreq, nbBit] = wavread([pathName, fileName]);
    if size(newSig, 2) > 1
      % multichannel wave, we keep only the first channel
      newSig = newSig(:, 1);
    end
    sig.sampFreq = sampFreq;
    sig.nbBit = nbBit;
    sig.ori = newSig;
    sig.mod = newSig;
    sig.real = true;
    
    resetMulaclab();
  end

  function saveModSig(objId, eventData)
    [fileName, pathName] = uiputfile('*.wav', 'Save modified signal');
    if fileName == 0
      % user pressed cancel, stop here
      return;
    end
    wavwrite(sig.mod, sig.sampFreq, sig.nbBit, [pathName, fileName]);
  end

  function openOriDecompo(objId, eventData)
    [fileName, pathName] = uigetfile('*.mat', 'Load frame parameters');
    if fileName == 0
      % user pressed cancel, stop here
      return;
    end
    data = load([pathName, fileName]);
    if ~(isfield(data, 'frame') && isfield(data, 'savedSig'))
      errordlg([fileName ' is not a valid signal decomposition file']);
      return;
    end
    
    frame = data.frame;
    sig = data.savedSig;
    sig.mod = sig.ori;
    
    resetMulaclab();
  end

  function saveOriDecompo(objId, eventData)
    [fileName, pathName] = uiputfile('*.mat', 'Export selection');
    if fileName == 0
      % user pressed cancel, stop here
      return;
    end
    savedSig = sig;
    rmfield(savedSig, 'mod');
    save([pathName, fileName], 'frame', 'savedSig');
  end

  function saveModDecompo(objId, eventData)
    [fileName, pathName] = uiputfile('*.mat', 'Export selection');
    if fileName == 0
      % user pressed cancel, stop here
      return;
    end
    savedSig = sig;
    savedSig.ori = savedSig.mod;
    rmfield(savedSig, 'mod');
    save([pathName, fileName], 'frame', 'savedSig');
  end

  function exportSel(objId, eventData)
    [fileName, pathName] = uiputfile('*.mat', 'Export selection');
    if fileName == 0
      % user pressed cancel, stop here
      return;
    end
    save([pathName, fileName], 'sel');
  end

  function exportSymbol(objId, eventData)
    [fileName, pathName] = uiputfile('*.png', 'Export symbol as...');
    if fileName == 0
      % user pressed cancel, stop here
      return;
    end
    cm = get(gui.mainFigId,'Colormap');
    oriInd = find(strcmp('originalMain', {visu.label}), 1);
    Clim = get(visu(oriInd).axesId,'CLim');
    oriIdx = round((flipud(coeff.oriC)-Clim(1))/(Clim(2)-Clim(1))*255)+1;
    oriExt = zeros([size(coeff.oriC),3]);
    oriExt(:,:,1) = reshape(cm(oriIdx,1),size(coeff.oriC));
    oriExt(:,:,2) = reshape(cm(oriIdx,2),size(coeff.oriC));
    oriExt(:,:,3) = reshape(cm(oriIdx,3),size(coeff.oriC));
    
% exporting imported symbols is not supported (does not make sence)   
%     symb = symbol.data(symbol.curSymb).val;
%     if(isempty(symb))
%        absSymb = flipud(abs(convSelToSymb()));
%     else
%        absSymb = flipud(abs(symb)); 
%     end

    absSymb = flipud(abs(convSelToSymb()));
    imwrite(oriExt,[pathName, fileName],'png','Alpha',absSymb,'bitdepth',16);
  end

  function importSel(objId, eventData)
    [fileName, pathName] = uigetfile('*.mat', 'Import selection');
    if fileName == 0
      % user pressed cancel, stop here
      return;
    end
    data = load([pathName, fileName]);
    if ~isfield(data, 'sel')
      errordlg([fileName ' is not a valid signal selection file']);
      return;
    end
    
    % as the saved selection might have been created using a longer signal,
    % or a signal with a higher sampling frequency, it could be wider that
    % the part of the time-frequency plane corresponding to the signal. 
    % We need to force the selection to stay in the right region
    fullPoly = fullSigPoly;
    for ind = 1:length(data.sel.lay)
      % intersection of the polygon with a polygon covering the whole
      % signal
      data.sel.lay(ind).poly = PolygonClip(data.sel.lay(ind).poly,...
        fullPoly, 1);
    end
    
    delSel;
    sel.lay(end+1:end+length(data.sel.lay)) = data.sel.lay;
    updateLayerList;
    set(gui.layListId, 'Value', 1);
    changeSelLay;
    drawAllSel;
  end

 function importSymbol(objId, eventData)
    [fileName, pathName] = uigetfile('*.png', 'Import symbol');
    if fileName == 0
      % user pressed cancel, stop here
      return;
    end
    [imag, ~, alpha] = imread([ pathName,fileName]);
    if(size(imag,3)>1)
      if isempty(alpha)
        errordlg(sprintf('%s does not contain alpha channel.',fileName));
        return;
      end
      mask = double(alpha)/double(max(alpha(:)));
    else
      mask = double(imag)/double(max(imag(:)));
    end

    if any(size(coeff.oriC)~=size(mask))
      errordlg([fileName ' is not a valid symbol image. The dimensions are different.']);
      return;
    end
    
    %limit
    mmss = 10^(default.zerotodb/20);
    mask(abs(mask)<mmss) = mmss;

    addSymbolListItem(flipud(mask));
    symbol.curSymb = numel(symbol.data);
    set(gui.symbListId, 'Value',symbol.curSymb);
    updateSymbolList();
    handleSymb();
  end

  % Menu callback function. When frame parameters are changed.
  function changeFrameDef(objId, eventData)
    newFrame = changeFrameDialog(frame);
    if ~isempty(newFrame)
      frame = newFrame;
      resetMulaclab();
    end
  end

  % When frame type is changed
  function changeFrameType(objId, eventData)
    resetSymbol(); 
     
    supportedFramesIdx = get(objId,'UserData'); 
    newFrame = supportedFrames{supportedFramesIdx};
    %newFrame = changeGaborFrameDialog(frame);
    if ~isempty(newFrame)
      frame = newFrame;
      resetMulaclab(); 
    end
  end

  function newFrame = changeFrameDialog(oldFrame)
    newFrame = [];
    
    margin = [10, 10];
    textSize = [170, 20];
    editSize = [80, 20];
    spaceSize = [10, 10];
    
    nbParam = numel(fieldnames(oldFrame.def));
    dialogPosition = [200, 200, ...
      2*margin(1)+textSize(1)+editSize(1)+spaceSize(1),...
      2*margin(2)+(nbParam+1)*textSize(2)+nbParam*spaceSize(1)];
    
    dialogId = dialog(...
      'Name', [oldFrame.type,' frame parameters'],...
      'Position', dialogPosition);
    
    buttonWidth = 50;
    buttonX = round((dialogPosition(3) -...
      (2*buttonWidth + spaceSize(1)))/2);
    
    buttonOkId = uicontrol(...
      'Parent', dialogId,...
      'Style', 'pushbutton',...
      'String', 'OK',...
      'Callback', @callbackOk,...
      'Position', [buttonX, margin(2), buttonWidth, textSize(2)]);
    
    
    buttonCancelId = uicontrol(...
      'Parent', dialogId,...
      'Style', 'pushbutton',...
      'String', 'Cancel',...
      'Callback', @callbackCancel,...
      'Position', [buttonX+buttonWidth+spaceSize(1), margin(2),...
        buttonWidth, textSize(2)]);
     
    oldFramedef = oldFrame.def;
    fieldNames = fieldnames(oldFramedef);
    fieldNames = fieldNames(end:-1:1);
    objIds = zeros(size(fieldNames));
    
    for ind = 1:length(fieldNames)
       posText = [margin+[0, ind*(spaceSize(2)+textSize(2))], textSize];
       posEdit = [margin+[0,...
          ind*(spaceSize(2)+textSize(2))]+[textSize(1)+spaceSize(1), 0],...
          editSize];
       
       propLabel = fieldNames{ind};
       if(exist('mulaclabDict','var')&&isfield(mulaclabDict,propLabel))
          propLabel = getfield(mulaclabDict,propLabel);
       end
          
        uicontrol(...
       'Parent', dialogId,...
       'Style', 'text',...
       'String', propLabel,...
       'Position', posText);
   
      strVal = getfield(oldFramedef,fieldNames{ind});
      if(isnumeric(strVal))
         strVal = num2str(strVal);
      end
       objIds(ind) = uicontrol(...
       'Parent', dialogId,...
       'Style', 'edit',...
       'String', strVal,...
       'Position', posEdit);
    end
        
    uiwait(dialogId);
    
    function callbackOk(objId, eventData)
      newFrame = oldFrame;
      for ind2 = 1:length(fieldNames)
         val = get(objIds(ind2),'string');
          if isempty(val)
            errordlg(sprintf('%s parameter is invalid.',fieldNames{ind2}));
            return;
          end
          if(isnumeric(getfield(oldFramedef,fieldNames{ind2})))
             newFrame.def = setfield(newFrame.def,fieldNames{ind2},round(str2num(val)));
          else
             newFrame.def = setfield(newFrame.def,fieldNames{ind2},val);
          end
      end
      close(dialogId);
    end
    
    function callbackCancel(objId, eventData)
      newFrame = [];
      close(dialogId);
    end
    
  end

  function loadFrameParam(objId, eventData)
    [fileName, pathName] = uigetfile('*.mat', 'Load frame parameters');
    if fileName == 0
      % user pressed cancel, stop here
      return;
    end
    data = load([pathName, fileName]);
    if ~isfield(data, 'frame')
      errordlg([fileName ' is not a valid signal frame file']);
      return;
    end
    frame = data.frame;
    coeff.ori = calculateCoeff(sig.ori);
    coeff.mod = calculateCoeff(sig.mod);
    updateVisu;
  end

  function saveFrameParam(objId, eventData)
    [fileName, pathName] = uiputfile('*.mat', 'Save frame parameters');
    if fileName == 0
      % user pressed cancel, stop here
      return;
    end
    save([pathName, fileName], 'frame');
  end

  function pos = buttonPos(indHori, indVert, shift) 
    % compute the position of a button (useful to align buttons)
    % indHori: horizontal index of the button for which we want the 
    % position (first element with index 1 on the left)
    % indVert: vertical index of the button for which we want the 
    % position (first element with index 1 at the bottom)
    
    pos = [shift(1) + (indHori-1)*(gui.buttonWidth+gui.horiDist)-1,...
      shift(2) + (indVert-1)*(gui.buttonHeight+gui.vertDist)-1,...
      gui.buttonWidth,...
      gui.buttonHeight];
    
  end

  function subSize = subpanelSize(nbRow)     
    subSize = [(gui.panelWidth - 4),...
      (2*gui.margin(2) + nbRow*gui.buttonHeight +...
      (nbRow-1) * gui.vertDist + gui.fontSize -4)];
  end

  function [nbFreq] = nbPlottedFreq()
     nbFreq = size(coeff.oriC, 1);
%     if sig.real
%       nbFreq = floor(size(coeff.ori, 1)/2)+1;
%     else
%       nbFreq = size(coeff.ori, 1);
%     end
  end

  function [] = resize(objId, eventData)
    figPos = get(gui.mainFigId, 'Position');
    xVisu = gui.panelWidth / figPos(3);
    % we need the follwing test to avoid a bug during creation
    if isinf(xVisu) || isnan(xVisu) 
      xVisu = 0.5;
    end
    
    widthVisu = (figPos(3)-gui.panelWidth) / figPos(3);
    % we need the follwing test to avoid a bug during creation
    if ~(widthVisu > 0)
      widthVisu = 0.5;
    end
    
    curY = 1;
    for ind = 1:length(visu)
      if visu(ind).visible
        curY = curY - visu(ind).height;
        visu(ind).position = [xVisu curY widthVisu visu(ind).height];
        set(visu(ind).uipanelId, 'Position', visu(ind).position);
        
        if ishandle(gui.ResizeVisuButtonId(ind))
          set(gui.ResizeVisuButtonId(ind), 'Position',...
            [figPos(3)-gui.visuResizeButtonWidth,...
              round(curY*figPos(4)-gui.visuResizeButtonHeight/2),...
              gui.visuResizeButtonWidth, ...
              gui.visuResizeButtonHeight]);
        end
      end
    end
  end

  function updateUndoData()
    undoData{end+1} = sel;
    set(gui.undoMenuId, 'Enable', 'on');
    redoData = {};
    set(gui.redoMenuId, 'Enable', 'off');
  end

  function undo(objId, eventData)
    % PI: Currently the implementation of undo/redo is very rudimentary. 
    % Only undo/redo of actions modifying the polygons of the selection sel
    % are taken into account. For example, zoom actions, or change of layer 
    % parameters (type of conversion, color of line, ...) are not taken 
    % into account. A more general handling could be implemented.
    % The memory consumption for undo could also be improved. Currently for
    % each action modifying the polygons, the whole old variable sel is
    % stored in a cell of undoData. It's generally contains more
    % information than stricly needed to undo as only one layer is modified
    % at a time. If improved efficiency is needed, it could be taken into
    % account to save minimally needed information. The number of undo
    % levels is currently unlimited, this could be also changed if needed
    if ~isempty(undoData)
      delSel;
      redoData{end+1} = sel;
      set(gui.redoMenuId, 'Enable', 'on');
      sel = undoData{end};
      undoData = undoData(1:end-1);
      if isempty(undoData)
        set(gui.undoMenuId, 'Enable', 'off');
      end
      drawAllSel;
      set(gui.layListId, 'Value', sel.curLay);
      updateLayerList;
    end
  end

  function redo(objId, eventData)
    % PI: see the PI for undo function
    if ~isempty(redoData)
      delSel;
      undoData{end+1} = sel;
      set(gui.undoMenuId, 'Enable', 'on');
      sel = redoData{end};
      redoData = redoData(1:end-1);
      if isempty(redoData)
        set(gui.redoMenuId, 'Enable', 'off');
      end
      drawAllSel;
      set(gui.layListId, 'Value', sel.curLay);
      updateLayerList;
    end
  end

% ______________________ VISUALIZATION FUNCTIONS __________________________

  function updateVisu(resetVisu, changeCLim)
   
    if nargin == 0
      resetVisu = false;
      changeCLim = false;
    end
    
    delete(gui.symbImageId);
    gui.symbImageId = [];
    
    % keep the properties of the axes if we're not resetting the
    % visualization
    axesProp = struct;
    if ~resetVisu
      for ind = 1:length(visu)
        if ~isempty(visu(ind).axesId)
          axesProp(ind).xLim = get(visu(ind).axesId, 'XLim');
          axesProp(ind).yLim = get(visu(ind).axesId, 'YLim');
          axesProp(ind).cLim = get(visu(ind).axesId, 'CLim');
        end
      end
    end
    
    % update the visualization of the signal
    linkXId = [];
    linkYId = [];
    linkCLimId=[];
    for ind = 1:length(visu)
      delete(visu(ind).uipanelId);
      visu(ind).uipanelId = [];
      visu(ind).axesId = [];
      if visu(ind).visible
        visu(ind).uipanelId = uipanel(gui.mainFigId, ...
          'Units', 'normalized', ...
          'Position', visu(ind).position);
 
        if visu(ind).linkX
          linkXId = [linkXId ind];
        end
        if visu(ind).linkY
          linkYId = [linkYId ind];
        end
        if visu(ind).linkCLim
          linkCLimId = [linkCLimId ind];
        end
      end
    end
    
    function updateOneVisu(ind)
         if visu(ind).visible
           visu(ind).axesId = visu(ind).updateFcn(visu(ind));

           % plot for the visualisation of play position during playback
           lim = axis(visu(ind).axesId);
           % xPos = convIndToAxesX((player.position-1) / frame.def.hop);
            hop = numel(sig.ori)/size(coeff.oriC,2);
            xPos = convIndToAxesX((player.position-1) / hop);
           hold(visu(ind).axesId, 'on');
           player.positionPlotId(ind) = plot(visu(ind).axesId, [xPos xPos],...
             lim(3:4), 'w'); % !!! give possibility to change color ?
           set(player.positionPlotId(ind), 'Visible', 'off', 'LineWidth', 2);
           hold(visu(ind).axesId, 'off');

           % the following insures the possibility to zoom out
           % completely even if we are already zoomed in
           axes(visu(ind).axesId);
           zoom(gui.mainFigId, 'reset');
        end
    end
    
    oriInd = find(strcmp('originalMain', {visu.label}), 1);
    updateOneVisu(oriInd);
    
    if length(linkCLimId) > 1
      link.CLim = get(visu(oriInd).axesId,'CLim');
    end
    range = 1:length(visu);
    range(oriInd) = [];
    for ind = range
       updateOneVisu(ind);
    end
     
    if length(linkXId) > 1
      link.XLim = linkprop(arrayfun(@(vEl) vEl.axesId,visu(linkXId)),'XLim');
    end
    if length(linkYId) > 1
      link.YLim = linkprop(arrayfun(@(vEl) vEl.axesId,visu(linkYId)),'YLim');
    end
    
    
    if ~resetVisu
      % we go in reverse order to update axes limits so that
      % all the linked axes come back to the limits of the original
      % spectrogram, even if some new linked axes are added
      for ind = length(visu):-1:1
        if visu(ind).visible
          try
            set(visu(ind).axesId, 'XLim', axesProp(ind).xLim);
          end
          try
            set(visu(ind).axesId, 'YLim', axesProp(ind).yLim);
          end
          if ~changeCLim
            try
              set(visu(ind).axesId, 'CLim', axesProp(ind).cLim);
            end
          end
        end
      end
    end
    
    drawAllSel;
    
    overviewInd = find(strcmp('originalOverview', {visu.label}), 1);
    if ~isempty(overviewInd)
      if visu(overviewInd).visible
        exploreOverview(visu(overviewInd).axesId);
      end
    end

    if ~resetVisu
      % restore the current tool
      changeTool([], [], gui.curTool);
    end
    
    drawResizeVisuButton;
    resize;    
  end

  function [axesId,C] = updateVisuCommon(coef,lim,visuData) 
    axesId = axes('Parent',visuData.uipanelId);
    if(lim)
       C = plotCoeff(coef,axesId,'clim',link.CLim);
    else
       C = plotCoeff(coef,axesId,'dynrange',visuData.dynamic);
    end   
    colorbar off;
  end

  function [axesId] = updateVisuOriSpec(visuData)
    [axesId, coeff.oriC] = updateVisuCommon(coeff.ori,false,visuData);
    title('Original signal');
  end

  function [axesId] = updateVisuModSpec(visuData)
   [axesId, coeff.modC] = updateVisuCommon(coeff.mod,visuData.linkCLim,visuData); 
%     if(visuData.linkCLim && ~strcmpi('originalOverview', {visu.label}))
%        plotCoeff(coeff.mod,axesId,'clim',link.CLim);
%     else
%        plotCoeff(coeff.mod,axesId,'dynrange',visuData.dynamic);
%     end
    %plotCoeff(coeff.mod,axesId,visuData.dynamic);
    % update the spectrogram of the original signal
%     spec = abs(coeff.mod(1:nbPlottedFreq,:)).^2; % spectrogram coefficients
%     maxSpec = max(spec(:));
%     imagesc(...
%       convIndToAxesX([1,size(spec,2)]),...
%       convIndToAxesY([1,size(spec,1)]),...
%       10*log10(spec+eps),...
%       'Parent', axesId, ...
%       10*log10(maxSpec)+[-visuData.dynamic 0]);
%     axis(axesId, 'xy');
%     xlabel('Time (s)');
%     ylabel('Frequency (Hz)');
    title('Modified signal');
  end

  function [axesId] = updateVisuOriOverview(visuData)
   axesId = updateVisuCommon(coeff.ori,visuData.linkCLim,visuData); 
   %axesId = updateVisuOriSpec(visuData);  
   % axesId = axes('parent',visuData.uipanelId);
   % plotCoeff(coeff.ori,axesId,visuData.dynamic);
    % update the spectrogram of the original signal
%     spec = abs(coeff.ori(1:nbPlottedFreq,:)).^2; %  spectrogram coefficients
%     maxSpec = max(spec(:));
%     imageId = imagesc(...
%       convIndToAxesX([1,size(spec,2)]),...
%       convIndToAxesY([1,size(spec,1)]),...
%       10*log10(spec+eps),...
%       'Parent', axesId, ...
%       10*log10(maxSpec)+[-visuData.dynamic 0]);
%     axis(axesId, 'xy');
%     xlabel('Time (s)');
%     ylabel('Frequency (Hz)');
    title('Overview of original signal');
  end

  function drawResizeVisuButton()
    for ind = 1:length(gui.ResizeVisuButtonId)
      if ishandle(gui.ResizeVisuButtonId(ind))
        delete(gui.ResizeVisuButtonId(ind));
      end
    end
    
    nbVisibleVisu = length(find([visu.visible]));
    nbButton = 0;
    gui.ResizeVisuButtonId = NaN(length(visu), 1);
    if nbVisibleVisu > 1
      resizeButtonIcon = imread([iconpath,'resizebutton.png']);
      for ind = 1:length(visu)
        if visu(ind).visible
          gui.ResizeVisuButtonId(ind) = uicontrol(...
            'Parent', gui.mainFigId,...
            'Style', 'pushbutton',...
            'String', '',...
            'TooltipString', 'Resize visualization',...
            'CData', resizeButtonIcon,...
            'BackgroundColor', gui.buttonBackgroundColor,...
            'Enable', 'off',...
            'buttonDownFcn', {@resizeVisu, ind});
          nbButton = nbButton + 1;
          if nbButton == nbVisibleVisu-1
            break;
          end
        end
      end
    end
  end

  function resizeVisu(objId, eventData, indVisu)
    
    initPos = get(gui.mainFigId, 'CurrentPoint');
    
    backupButtonMotionFcn = get(gui.mainFigId, 'WindowButtonMotionFcn');
    backupButtonUpFcn = get(gui.mainFigId, 'WindowButtonUpFcn');
    
    figPos = get(gui.mainFigId, 'Position');
    
    buttonPos = get(gui.ResizeVisuButtonId(indVisu), 'Position');
    
    % find the index of the visible visualization next to the one with 
    % index indVisu
    indNextVisu = 0;
    for ind = indVisu+1:length(visu)
      if visu(ind).visible
        indNextVisu = ind;
        break;
      end
    end
    
    initHeight = [visu(indVisu).height, visu(indNextVisu).height];
    
    % As the zoom and pan mode of matlab doesn't allow to change the
    % WindowButtonUpFcn, I temporarly change the tool
    curTool = gui.curTool;
    changeTool(objId, eventData, 'freehand');
    
    set(gui.mainFigId,...
      'WindowButtonMotionFcn', @resizeVisuButtonMotionFcn);
    set(gui.mainFigId,...
      'WindowButtonUpFcn', @resizeVisuButtonUpFcn);
    
    function resizeVisuButtonMotionFcn(objId, eventData)
      pos = get(gui.mainFigId, 'CurrentPoint');
      shift = (pos(2) - initPos(2)) / figPos(4);
      if (initHeight(1) - shift) > 0.01 && (initHeight(2) + shift) > 0.01
        visu(indVisu).height = initHeight(1) - shift;
        visu(indNextVisu).height = initHeight(2) + shift;
      end
      resize;
    end
    
    function resizeVisuButtonUpFcn(objId, eventData)
      set(gui.mainFigId, 'WindowButtonMotionFcn', backupButtonMotionFcn);
      set(gui.mainFigId, 'WindowButtonUpFcn', backupButtonUpFcn);
      % restore the tool
      changeTool(objId, eventData, curTool);
    end
  end

  function toggleModVisu(objId, eventData)
    oriInd = find(strcmp('originalMain', {visu.label}), 1);
    modInd = find(strcmp('modified', {visu.label}), 1);
    
    if visu(modInd).visible
      visu(modInd).visible = false;
      set(gui.visuMenu.showModId, 'Checked', 'off');
      visu(oriInd).height = visu(oriInd).height + visu(modInd).height;
    else
      visu(modInd).visible = true;
      set(gui.visuMenu.showModId, 'Checked', 'on');
      visu(modInd).height = visu(oriInd).height / 2;
      visu(oriInd).height = visu(oriInd).height / 2;
    end
    updateVisu;
    resize;
  end

  function toggleOverviewVisu(objId, eventData)
    oriInd = find(strcmp('originalMain', {visu.label}), 1);
    overviewInd = find(strcmp('originalOverview', {visu.label}), 1);
    
    if visu(overviewInd).visible
      visu(overviewInd).visible = false;
      set(gui.visuMenu.showOverviewId, 'Checked', 'off');
      visu(oriInd).height = visu(oriInd).height + visu(overviewInd).height;
    else
      visu(overviewInd).visible = true;
      set(gui.visuMenu.showOverviewId, 'Checked', 'on');
      visu(overviewInd).height = visu(oriInd).height * 1/3;
      visu(oriInd).height = visu(oriInd).height * 2/3;
    end
    updateVisu;
    resize;
  end



% ______________________ AUDIOPLAYER FUNCTIONS ____________________________

  function initializePlayer()
    
    player.ori = audioplayer(sig.ori, sig.sampFreq);
    player.mod = audioplayer(sig.mod, sig.sampFreq);
    
    set(player.ori, 'StartFcn', @playerStartFcn);
    set(player.ori, 'StopFcn', @playerStopFcnOri);
    set(player.ori, 'TimerFcn', @updatePlayPosition);

    set(player.mod, 'StartFcn', @playerStartFcn);
    set(player.mod, 'StopFcn', @playerStopFcnMod);
    set(player.mod, 'TimerFcn', @updatePlayPosition);

    player.selected = 'ori';
    player.loop = false;
    player.position = 1;
    player.forceStop = false;
    
    if ishandle(gui.player.buttongroupId)
      set(gui.player.buttongroupId,...
        'SelectedObject', gui.player.buttonOriId);
    end
  end


  function drawAudioplayer(uipanelId)    
    subPanelHeight = 2*gui.textHeight + 2*gui.marginSub(2) +...
      gui.vertDist+gui.fontSize+2;
    
    gui.player.buttongroupId = uibuttongroup(...
      'Parent', uipanelId,...
      'Title', 'Signal',...
      'TitlePosition', 'centertop',...
      'FontSize', gui.fontSize, ...
      'Units', 'pixels',...
      'Position', [1, 1, gui.panelWidth-4, subPanelHeight],...
      'SelectionChangeFcn', @switchAudioplayer);
      
    gui.player.buttonOriId = uicontrol(...
      'Parent', gui.player.buttongroupId,...
      'Style', 'radiobutton',...
      'FontSize', gui.fontSize, ...
      'String', ' Original',...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'Tag', 'ori',...
      'Position', [gui.marginSub(1),...
        gui.textHeight+gui.marginSub(2)+gui.vertDist,...
        gui.panelWidth-4-2*gui.marginSub(2),...
        gui.textHeight]);
    
    gui.player.buttonModId = uicontrol(...
      'Parent', gui.player.buttongroupId,...
      'Style', 'radiobutton',...
      'FontSize', gui.fontSize, ...
      'String', ' Modified',...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'Tag', 'mod',...
      'Position', [gui.marginSub,...
        gui.panelWidth-4-2*gui.marginSub(2),...
        gui.textHeight]);
    
    playIcon = imread([iconpath,'play.png']);
    
    yPos = subPanelHeight;
        
    buttonPlayPauseId = uicontrol(...
      'Parent', uipanelId,...
      'Style', 'pushbutton',...
      'String', '',...
      'TooltipString', 'Play/pause audio',...
      'CData', playIcon,...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'Position', buttonPos(1, 1, gui.margin + [0, yPos]),...
      'CallBack',@playPauseAudioplayer);

    stopIcon = imread([iconpath,'stop.png']);
    
    buttonStopId = uicontrol(...
      'Parent', uipanelId,...
      'Style', 'pushbutton',...
      'String', '',...
      'TooltipString', 'Stop audio',...
      'CData', stopIcon,...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'Position', buttonPos(2, 1, gui.margin + [0, yPos]),...
      'CallBack',@stopAudioplayer);
    
    loopIcon = imread([iconpath,'loop.png']);
    
    buttonLoopId = uicontrol(...
      'Parent', uipanelId, ...
      'Style', 'togglebutton',...
      'String', '',...
      'TooltipString', 'Loop audio',...
      'CData', loopIcon,...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'Position', buttonPos(3, 1, gui.margin + [0, yPos]),...
      'CallBack',@loopAudioplayer,...
      'Min', false,...
      'Max', true);
    
  end

  function switchAudioplayer(objId, eventData) 
    player.selected=get(get(objId, 'SelectedObject'),'Tag');
  end

  function playPauseAudioplayer(objId, eventData)
    if isplaying(player.ori)
      player.forceStop = true;
      pause(player.ori);
    elseif isplaying(player.mod)
      player.forceStop = true;
      pause(player.mod);
    else
      player.forceStop = false;
      switch player.selected
        case 'ori'
          play(player.ori, player.position);
        case 'mod'
          play(player.mod, player.position);
      end
    end
  end


  function stopAudioplayer(objId, eventData)
    player.forceStop = true;
    stop(player.ori);
    stop(player.mod);
    player.position = 1;
  end

  function playerStartFcn(objId, eventData)
    for ind = 1:length(player.positionPlotId)
      if visu(ind).visible
        set(player.positionPlotId(ind), 'Visible', 'on');
      end
    end
  end

  function playerStopFcnOri(objId, eventData)
    for ind = 1:length(player.positionPlotId)
      if visu(ind).visible
        set(player.positionPlotId(ind), 'Visible', 'off');
      end
    end   
    
    player.position = get(player.ori, 'CurrentSample');
    
    if player.loop && ~player.forceStop
      pause(0.1); %!!! needed to avoid java error
      play(player.ori, player.position);
    end
  end

  function playerStopFcnMod(objId, eventData)
    for ind = 1:length(player.positionPlotId)
      if visu(ind).visible
        set(player.positionPlotId(ind), 'Visible', 'off');
      end
    end
    
    player.position = get(player.mod, 'CurrentSample');
    
    if player.loop && ~player.forceStop
      pause(0.1); %!!! needed to avoid java error
      play(player.mod, player.position);
    end
  end

  function loopAudioplayer(objId, eventData)
    player.loop = get(objId, 'Value');
  end

  function updatePlayPosition(objId, eventData)
    % pos = convIndToAxesX((get(objId, 'CurrentSample')-1) / frame.def.hop);
    hop = numel(sig.ori)/size(coeff.oriC,2);
    pos = convIndToAxesX((get(objId, 'CurrentSample')-1) / hop);
    
    for ind = 1:length(player.positionPlotId)
      if visu(ind).visible
        set(player.positionPlotId(ind), 'XData', [pos pos]);
      end
    end
  end

  function updatePlayerMod()
    player.mod = audioplayer(sig.mod, sig.sampFreq);
    set(player.mod, 'StartFcn', @playerStartFcn);
    set(player.mod, 'StopFcn', @playerStopFcnMod);
    set(player.mod, 'TimerFcn', @updatePlayPosition);
  end

% ____________________ SELECTION TOOLS FUNCTIONS __________________________

  function drawSelectionTool(uipanelId)
    subSize = subpanelSize(1);
    
    processingToolPanelId = uipanel(...
      'Parent', uipanelId,...
      'Title', 'Apply',...
      'TitlePosition', 'centertop',...
      'FontSize', gui.fontSize, ...
      'Units', 'pixels',...
      'Position', [1, 1, subSize]);
    
    % draw processing tools subpanel
    drawProcessingTool(processingToolPanelId);
    
    heightLayerPanel = 115;
    
    layerPanelId = uipanel(...
      'Parent', uipanelId,...
      'Title', 'Layers',...
      'TitlePosition', 'centertop',...
      'FontSize', gui.fontSize, ...
      'Units', 'pixels',...
      'Position', [1, 1+subSize(2), subSize(1), heightLayerPanel]);
    
    % draw layer panel
    drawLayerPanel(layerPanelId);
    
    buttongroupId = uibuttongroup(...
      'Parent', uipanelId,...
      'Title', 'Mode',...
      'TitlePosition', 'centertop',...
      'FontSize', gui.fontSize, ...
      'Units', 'pixels',...
      'Position', [1, subSize(2)+heightLayerPanel+1, subSize],...
      'SelectionChangeFcn', @switchSelMode);
    
    unionIcon = imread([iconpath,'union.png']);
    
    buttonUnionId = uicontrol(...
      'Parent',buttongroupId,...
      'Style', 'togglebutton',...
      'String', '',...
      'TooltipString', 'Union',...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'CData', unionIcon,...
      'Tag', 'union',...
      'Position', buttonPos(1, 1, gui.marginSub));
    
    intersectionIcon = imread([iconpath,'intersection.png']);
    
    buttonInterId = uicontrol(...
      'Parent',buttongroupId,...
      'Style', 'togglebutton',...
      'String', '',...
      'TooltipString', 'Intersection',...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'CData', intersectionIcon,...
      'Tag', 'intersection',...
      'Position', buttonPos(2, 1, gui.marginSub));
    
    differenceIcon = imread([iconpath,'difference.png']);
    
    buttonDiffId = uicontrol(...
      'Parent',buttongroupId,...
      'Style', 'togglebutton',...
      'String', '',...
      'TooltipString', 'Set difference',...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'CData', differenceIcon,...
      'Tag', 'difference',...
      'Position', buttonPos(3, 1, gui.marginSub));
    
    % reserved space to draw parameters of the tools
    gui.toolPanelId = uipanel(...
      'Parent', uipanelId,...
      'Title', 'Tools',...
      'TitlePosition', 'centertop',...
      'FontSize', gui.fontSize,...
      'Units', 'pixels',...
      'Position', [1, 2*subSize(2)+heightLayerPanel+1,...
        subSize(1),...
        2*gui.marginSub(2) + gui.vertDist + gui.textHeight + ...
        gui.buttonHeight + gui.fontSize+2]);
    
      
      
    
    freeHandIcon = imread([iconpath,'freehand.png']);
    
    buttonFreehandId = uicontrol(...
      'Parent', gui.toolPanelId,...
      'HandleVisibility', 'off',...
      'Style', 'togglebutton',...
      'String', '',...
      'TooltipString', 'Free-hand selection',...
      'CData', freeHandIcon,...
      'Min', false,...
      'Max', true,...
      'Value', false,...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'Position', buttonPos(1, 1,...
        gui.marginSub+[0,...
        gui.marginSub(2)+gui.textHeight]),...
      'CallBack', {@changeTool, 'freehand'});
    
    gui.tool(end+1).buttonId = buttonFreehandId;
    gui.tool(end).name = 'freehand';
    gui.tool(end).function = @selecFreehand;
    
    levelIcon = imread([iconpath,'magicwand.png']);
    
    buttonLevelId = uicontrol(...
      'Parent', gui.toolPanelId,...
      'HandleVisibility', 'off',...
      'Style', 'togglebutton',...
      'String', '',...
      'TooltipString', 'Magic wand (Level selection)',...
      'CData', levelIcon,...
      'Min', false,...
      'Max', true,...
      'Value', false,...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'Position', buttonPos(2, 1,...
        gui.marginSub+[0,...
        gui.marginSub(2)+gui.textHeight]),...
      'CallBack', {@changeTool, 'level'});
    
    gui.tool(end+1).buttonId = buttonLevelId;
    gui.tool(end).name = 'level';
    gui.tool(end).function = @selecLevel;
    gui.tool(end).param.name = 'Tolerance';
    gui.tool(end).param.val = '10'; % value of the tolerance in dB
    
    
    freeHandIcon = imread([iconpath,'subbandsel.png']);
    
    buttonSubbandId = uicontrol(...
      'Parent', gui.toolPanelId,...
      'HandleVisibility', 'off',...
      'Style', 'togglebutton',...
      'String', '',...
      'TooltipString', 'Subband selection',...
      'CData', freeHandIcon,...
      'Min', false,...
      'Max', true,...
      'Value', false,...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'Position', buttonPos(3, 1,...
        gui.marginSub+[0,...
        gui.marginSub(2)+gui.textHeight]),...
      'CallBack', {@changeTool, 'subband'});
    
    gui.tool(end+1).buttonId = buttonSubbandId;
    gui.tool(end).name = 'subband';
    gui.tool(end).function = @selecSubband;
  end

  function changeTool(objId, eventData, toolName)
    if getMatlabVersion < 7.3
      % As we cannot force automatic update of rectangle in overview when
      % zooming and paning, we will at least force it when we change tool
      updateExploreRect;
    end
    zoom(gui.mainFigId, 'off');
    pan(gui.mainFigId, 'off');    
    for ind = 1:length(gui.tool)
      set(gui.tool(ind).buttonId, 'Value', false);
    end
    
    toolNameList = {gui.tool.name}; % cell array listing the tool names
    curInd = find(strcmp(toolName, toolNameList));
    set(gui.tool(curInd).buttonId, 'Value', true);
    gui.curTool = gui.tool(curInd).name;    
    gui.tool(curInd).function(curInd);
    drawToolParam(curInd);
  end

  function drawToolParam(toolInd)    
    
    child = findobj(gui.toolPanelId);
    child = setdiff(child, gui.toolPanelId);
    delete(child);    
    
    for ind = 1:length(gui.tool(toolInd).param)
      uicontrol(...
        'Parent', gui.toolPanelId,...
        'FontSize', gui.fontSize, ...
        'Style', 'text',...
        'String', gui.tool(toolInd).param(ind).name,...
        'Position', [gui.marginSub(1),...
          gui.marginSub(2)+(ind-1)*(gui.textHeight+gui.vertDist),...
          gui.textWidth, gui.textHeight]) ;
      uicontrol(...
        'Parent', gui.toolPanelId,...
        'FontSize', gui.fontSize, ...
        'Style', 'edit',...
        'BackgroundColor', gui.buttonBackgroundColor,...
        'String', gui.tool(toolInd).param(ind).val,...
        'Position', [gui.marginSub(1)+gui.textWidth,...
          gui.marginSub(2)+(ind-1)*(gui.textHeight+gui.vertDist),...
          gui.editWidth, gui.textHeight],...
        'Callback', {@changeToolParam, toolInd, ind});
    end
  end

  function changeToolParam(objId, eventData, toolInd, paramInd)
    gui.tool(toolInd).param(paramInd).val = get(objId, 'String');
  end

  function selecFreehand(toolInd)
    % toolInd is not used but needed as I pass a parameter to this 
    % function !!! precise this

    figId = gui.mainFigId;
    oriInd = find(strcmp('originalMain', {visu.label}), 1);
    axesId = visu(oriInd).axesId;
    imageId = findobj(axesId, 'Type', 'Image');
    
    set(imageId, 'ButtonDownFcn', @selecFreeHandButtonDown);
    
    poly.x = [];
    poly.y = [];

    
    init = true; % boolean to know if we are drawing the first point or not
    
    lineId = [];
    
    function selecFreeHandButtonDown(objId, eventData)
      if strcmp(get(figId,'SelectionType'),'normal')
        if init
          point = get(axesId,'CurrentPoint');
          %disp(sprintf('[%d,%d]\n',point(1,1),point(1,2)));
          poly.x = [point(1,1)];
          poly.y = [point(1,2)];
          poly.hole = false;
          lineId = line(...
            'XData', poly.x,...
            'YData', poly.y,...
            'LineWidth', sel.lay(sel.curLay).lineWidth,...
            'Color', sel.lay(sel.curLay).color,...
            'Marker', sel.lay(sel.curLay).marker,...
            'MarkerSize', sel.lay(sel.curLay).markerSize,...
            'LineStyle', sel.lay(sel.curLay).lineStyle);
          init = false;
          set(figId,'WindowButtonDownFcn',@selecFreeHandButtonDown);
          set(figId,'WindowButtonMotionFcn',@selecFreeHandMotionDown);
          set(figId,'WindowButtonUpFcn',@selecFreeHandButtonUp);
        else
          point = get(axesId,'CurrentPoint');
          %disp(sprintf('[%d,%d]\n',point(1,1),point(1,2)));
          poly.x = [poly.x; point(1,1)];
          poly.y = [poly.y; point(1,2)];
          set(lineId, 'XData',poly.x , 'YData', poly.y);
          set(figId,'WindowButtonMotionFcn',@selecFreeHandMotionDown);
          drawnow;
        end
      end
    end

    function selecFreeHandMotionDown(objId, eventData)
      point = get(axesId,'CurrentPoint');
      % disp(sprintf('[%d,%d]\n',point(1,1),point(1,2)));
      poly.x = [poly.x; point(1,1)];
      poly.y = [poly.y; point(1,2)];
      set(lineId, 'XData', poly.x, 'YData', poly.y);
      drawnow;
    end

    function selecFreeHandMotionUp(objId, eventData)
      point = get(axesId,'CurrentPoint');
      xData = [poly.x; point(1,1)];
      yData = [poly.y; point(1,2)];
      set(lineId, 'XData', xData, 'YData', yData);
      drawnow;
    end

    function selecFreeHandButtonUp(objId, eventData)
      switch get(figId,'SelectionType')
        case 'normal'
          set(figId,'WindowButtonMotionFcn',@selecFreeHandMotionUp);
        case {'open', 'alt'}
          set(figId,'WindowButtonMotionFcn','');
          set(figId,'WindowButtonDownFcn','');
          set(figId,'WindowButtonUpFcn','');
          init = true;
          
          delete(lineId);
          
          % !!! we could simplify the polygon here to have fastest
          % computation 
          
          combineSel(poly);
          drawCurSel;
      end
    end
  end

  function selecLevel(toolInd)

    oriInd = find(strcmp('originalMain', {visu.label}), 1);
    imageId = findobj(visu(oriInd).axesId, 'Type', 'Image');
    
    set(imageId, 'ButtonDownFcn', @selecLevelButtonDown);
    
    function selecLevelButtonDown(objId, eventData)      
      point = get(visu(oriInd).axesId,'CurrentPoint');

      indTime = round(convAxesToIndX(point(1,1)));
      indFreq = ceil(convAxesToIndY(point(1,2)));
      %disp(sprintf('%i, %i',point(1,2),indFreq));

      % !!! sometime log of 0, do something here
      %spec = 20*log10(abs(coeff.ori(1:nbPlottedFreq,:))+eps);

      refLevel = coeff.oriC(indFreq, indTime);

      levelTol = str2num(gui.tool(toolInd).param.val);
      if isempty(levelTol) || levelTol < 0
        errordlg(['Tolerance parameter for magic wand is invalid, '...
          'it must be a positive number in dB']);
        return;
      end
      
      lowLevel = refLevel - levelTol;
      highLevel = refLevel + levelTol;

      mask = (coeff.oriC >= lowLevel) & (coeff.oriC <= highLevel);
      mask = bwselect(mask, indTime, indFreq, 4);

      poly = mask2poly(mask, indTime, indFreq);
      
      % !!! see if we simplify the selection or not
      
      for ind = 1:length(poly)
        poly(ind).x = convIndToAxesX(floor(poly(ind).x));
        poly(ind).y = convIndToAxesY(floor(poly(ind).y));
      end

      combineSel(poly);
      drawCurSel;
    end
    
  end

  function selecSubband(toolInd)
    oldShowSymb = [];


    oriInd = find(strcmp('originalMain', {visu.label}), 1);
    imageId = findobj(visu(oriInd).axesId, 'Type', 'Image');
    figId = gui.mainFigId;
    
    set(imageId, 'ButtonDownFcn', @selecSubbandButtonDown);
    
    selectedInFreq = zeros(size(coeff.oriC,1));
    startInFreq = 0;
    function selecSubbandButtonDown(objId, eventData)
      oldShowSymb =  gui.showSymb;
      gui.showSymb = false;
      handleSymb();
      point = get(visu(oriInd).axesId,'CurrentPoint');
      startInFreq = ceil(convAxesToIndY(point(1,2)));
      selecSubbandMotionDown(objId, eventData);
      set(figId, 'WindowButtonMotionFcn', @selecSubbandMotionDown);
      set(figId,'WindowButtonUpFcn',@selecSubbandButtonUp);
    end
    
    function selecSubbandMotionDown(objId, eventData)
      point = get(visu(oriInd).axesId,'CurrentPoint');
      indFreq = ceil(convAxesToIndY(point(1,2)));
       if(indFreq>size(coeff.oriC,1) || indFreq<1 )
         return;
      end

      selDif = indFreq - startInFreq;
      if selDif>0
         range = startInFreq:indFreq;
      else
         range = indFreq:startInFreq;
         range = range(end:-1:1);
      end
      if(all(selectedInFreq(range)==1))
      return;
      end
      

         poly.x = [1; size(coeff.oriC,2); size(coeff.oriC,2); 1];
         poly.y = [range(1)+0.5;range(1)+0.5;range(end)-0.5;range(end)-0.5];
         poly.hole = false;
      
         for ind = 1:length(poly)
           poly(ind).x = convIndToAxesX(floor(poly(ind).x));
           poly(ind).y = convIndToAxesY(floor(poly(ind).y));
         end

      combineSel(poly);
      drawCurSel;
      selectedInFreq(range) = 1;
    end
    
    function selecSubbandButtonUp(objId, eventData)
       gui.showSymb = oldShowSymb;
       handleSymb();
       set(figId,'WindowButtonMotionFcn','');
       set(figId,'WindowButtonUpFcn','');
       selectedInFreq = zeros(size(coeff.oriC,1));
    end
    
    set(figId,'WindowButtonMotionFcn','');
    set(figId,'WindowButtonUpFcn','');

  end

  function [poly] = mask2poly(mask, indTime, indFreq)
    % Convert region mask to region of interest (ROI) polygon
    % !!! I'm not shure that this algorithm works in every case

    maskRef = mask;
    mask = zeros(size(mask, 1)+2, size(mask, 2)+2);
    mask(2:end-1,2:end-1) = maskRef;

    [ind1, ind2] = find(mask);

    % position of segments on the edge for first dimension
    seg1 = sparse(size(mask,1)+1, size(mask,2)); 
    % position of segments on the edge for second dimension
    seg2 = sparse(size(mask,1), size(mask,2)+1) ;

    for n = 1:length(ind1)
      if ~mask(ind1(n)-1, ind2(n))
        seg1(ind1(n), ind2(n)) = 1;
      end
      if ~mask(ind1(n)+1, ind2(n))
        seg1(ind1(n)+1, ind2(n)) = 1;
      end
      if ~mask(ind1(n), ind2(n)-1)
        seg2(ind1(n), ind2(n)) = 1;
      end
      if ~mask(ind1(n), ind2(n)+1)
        seg2(ind1(n), ind2(n)+1) = 1;
      end
    end

    ind = 0;
    
    while nnz(seg1)>1
      curve = chainSeg;
      if ~isempty(curve)
        ind = ind+1;
        poly(ind).x = curve(:,2);
        poly(ind).y = curve(:,1);
        poly(ind).hole = true;
      end
    end
    
    poly(1).hole = false;
    
    function curve = chainSeg()
    % chaining of segments

      % choose one start segment
      [ind1, ind2] = find(seg1);
      ind1 = ind1(1);
      ind2 = ind2(1);
      curve = [ind1, ind2; ind1, ind2+1];

      seg1(ind1, ind2) = 0;

      % pos is a variable to remember in which direction we are
      % progressing
      % [1 1] if we're going in a growing index number
      % [-1 0] otherwise
      pos = [1, 1];

      % precise if the last added segment comes from seg1 or seg2
      last1 = true; 

      while true
        if last1
          % the last segment was from dimension 1
          if seg1(ind1, ind2+pos(1))
            ind2 = ind2+pos(1);
            seg1(ind1, ind2) = 0;
            curve = [curve; ind1, ind2+pos(2)];
          elseif seg2(ind1-1, ind2+pos(2))
            ind1 = ind1-1;
            ind2 = ind2+pos(2);
            pos = [-1, 0];
            seg2(ind1, ind2) = 0;
            curve = [curve; ind1, ind2];
            last1 = false;
          elseif seg2(ind1, ind2+pos(2))
            ind2 = ind2+pos(2);
            pos = [1, 1];
            seg2(ind1, ind2) = 0;
            curve = [curve; ind1+1, ind2];
            last1 = false;
          else
            break;
          end
        else
          % the last segment was from dimension 2
          if seg2(ind1+pos(1), ind2)
            ind1 = ind1+pos(1);
            seg2(ind1, ind2) = 0;
            curve = [curve; ind1+pos(2), ind2];
          elseif seg1(ind1+pos(2), ind2-1)
            ind1 = ind1+pos(2);
            ind2 = ind2-1;
            pos = [-1, 0];
            seg1(ind1, ind2) = 0;
            curve = [curve; ind1, ind2];
            last1 = true;
          elseif seg1(ind1+pos(2), ind2)
            ind1 = ind1+pos(2);
            pos = [1, 1];
            seg1(ind1, ind2) = 0;
            curve = [curve; ind1, ind2+1];
            last1 = true;
          else
            break;
          end
        end
      end
      
      curve = curve-1.5;
      
      if size(curve, 1)==2
        % the first segment couldn't be linked to some other segment,
        % there is no polygon
        curve = [];
      end
      
    end
    
  end

  function combineSel(poly)
    if ~isempty(sel.lay(sel.curLay).poly)
      % there is already a selection that should be combined with the new
      % one
      
      switch sel.mode
        case 'union'
          newPoly = PolygonClip(sel.lay(sel.curLay).poly, poly, 3);
        case 'intersection'
          newPoly = PolygonClip(sel.lay(sel.curLay).poly, poly, 1);
        case 'difference'
          newPoly = PolygonClip(sel.lay(sel.curLay).poly, poly, 0);
      end
      
      % we must remove the selection plot before replacing the selection
      % data
      delCurPoly;
      sel.lay(sel.curLay).poly = newPoly;   
    else
      % the selection was empty
      if strcmp(sel.mode, 'union')
        updateUndoData;
        sel.lay(sel.curLay).poly = poly;
      end
    end
  end

  function delCurPoly()
    updateUndoData;
    if ~isempty(sel.lay(sel.curLay).poly)
      for ind = 1:length(sel.lay(sel.curLay).poly)
        delete(sel.lay(sel.curLay).poly(ind).id);
      end
      sel.lay(sel.curLay).poly = [];
    end
  end

  function delAllPoly()
    updateUndoData;
    for indLay = 1:length(sel.lay)
      if ~isempty(sel.lay(indLay).poly)
        for ind = 1:length(sel.lay(indLay).poly)
          delete(sel.lay(indLay).poly(ind).id);
        end
        sel.lay(indLay).poly = [];
      end
    end
  end

  function delSel() % delete graphical selection (but keep the polygons)
    for indLay = 1:length(sel.lay)
      if ~isempty(sel.lay(indLay).poly)
        for ind = 1:length(sel.lay(indLay).poly)
          delete(sel.lay(indLay).poly(ind).id);
        end
      end
    end
  end

  function resetSel() % erase sel and put it to default value
    sel = default.sel;
    if isfield(gui, 'layListId')
      if ishandle(gui.layListId)
        set(gui.layListId, 'Value', sel.curLay);
        drawLayParam(sel.curLay);
        updateLayerList;
      end
    end
  end

  function drawCurSel()
    handleSymb;
    oriInd = find(strcmp('originalMain', {visu.label}), 1);
    hold(visu(oriInd).axesId, 'on');
    if ~isempty(sel.lay(sel.curLay).poly)
      for ind = 1:length(sel.lay(sel.curLay).poly)
        sel.lay(sel.curLay).poly(ind).id = line(...
          'XData', [sel.lay(sel.curLay).poly(ind).x;...
            sel.lay(sel.curLay).poly(ind).x(1)],...
          'YData', [sel.lay(sel.curLay).poly(ind).y;...
            sel.lay(sel.curLay).poly(ind).y(1)],...
          'LineWidth', sel.lay(sel.curLay).lineWidth,...
          'Color', sel.lay(sel.curLay).color,...
          'Marker', sel.lay(sel.curLay).marker,...
          'MarkerSize', sel.lay(sel.curLay).markerSize,...
          'LineStyle', sel.lay(sel.curLay).lineStyle,...
          'Parent', visu(oriInd).axesId);
      end
    end
    hold(visu(oriInd).axesId, 'off');
  end

  function drawAllSel()
     
    handleSymb;
    oriInd = find(strcmp('originalMain', {visu.label}), 1);
    hold(visu(oriInd).axesId, 'on');
    for indLay = 1:length(sel.lay)
      if ~isempty(sel.lay(indLay).poly)
        for ind = 1:length(sel.lay(indLay).poly)
          sel.lay(indLay).poly(ind).id = line(...
            'XData', [sel.lay(indLay).poly(ind).x;...
              sel.lay(indLay).poly(ind).x(1)],...
            'YData', [sel.lay(indLay).poly(ind).y;...
              sel.lay(indLay).poly(ind).y(1)],...
            'LineWidth', sel.lay(indLay).lineWidth,...
            'Color', sel.lay(indLay).color,...
            'Marker', sel.lay(indLay).marker,...
            'MarkerSize', sel.lay(indLay).markerSize,...
            'LineStyle', sel.lay(indLay).lineStyle,...
            'Parent', visu(oriInd).axesId);
        end
      end
    end
    hold(visu(oriInd).axesId, 'off');
  end

  function handleSymb()
    delete(gui.symbImageId);
    gui.symbImageId = [];
    if gui.showSymb
      showSymb();
      %gui.showSymb = false;
    end
  end

  function showSymb()
    % plot symbol of the multiplier
    oriInd = find(strcmp('originalMain', {visu.label}), 1);
    axesId = visu(oriInd).axesId;
    hold(axesId, 'on');
    
    symb = symbol.data(symbol.curSymb).val;
    if(isempty(symb))
       absSymb = abs(convSelToSymb());
    else
%        if export.limitXaxesRes
%         absSymb=interp1(linspace(0,1,size(symb,2)),abs(symb).',linspace(0,1,export.xLim),'nearest');
%         absSymb=absSymb.';
%        else
         absSymb = abs(symb); 
%       end
    end

     if(symbol.data(symbol.curSymb).invert)
       absSymb = max([1,max(absSymb(:))])-absSymb;
    end
    
    % Symbol is represented using transparency, with two different
    % colors for absolute values of the symbol in [0,1] and for values > 1. 
    % Transparency is total for absolute value of 1, and opacity
    % grows when going away from 1
    
    % !!! these colors should be set as variables that the user can set
    color1 = [0 0 0]; % black
    color2 = [1 1 1]; % white
    
    % !!! I should probably use a dB scale for transparency
    
    nbTime = size(coeff.oriC, 2);
    % construct the image with white for value > 1 and black in other parts
   
    tempIm = repmat(reshape(color1, [1 1 3]), [nbPlottedFreq nbTime 1]);
    [ind1, ind2] = find(absSymb(1:nbPlottedFreq, :) > 1);
    
    if ~isempty(ind1)
      for ind3 = 1:size(tempIm, 3)
        ind = sub2ind(size(tempIm), ind1, ind2, ind3 * ones(size(ind1)));
        tempIm(ind) = color2(ind3);
      end
    end

    
    % construct transparency data
    
    % initialisation using formula valid for values in [0, 1]
    tempAlpha = (1 - absSymb(1:nbPlottedFreq, :));
    
    % correction for values > 1 % !!! I don't really know what to use here
    tempAlpha(ind1, ind2) = 1 - 1./absSymb(ind1, ind2);
    
    gui.symbImageId = image(...
      convIndToAxesX([0, nbTime]),...
      convIndToAxesY([0.5, nbPlottedFreq-0.5]),...
      tempIm,...
      'alphaData',gui.symbOpacity * tempAlpha,...
      'Parent', axesId);
    hold(axesId, 'off');
    % restore the current tool
    changeTool([], [], gui.curTool);
  end
  
  function switchSelMode(objId, eventData)
    sel.mode = get(get(objId, 'SelectedObject'),'Tag');
  end

  function drawProcessingTool(uipanelId)
  % creation of apply button
  applyIcon = imread([iconpath,'apply.png']);
  
  applyId = uicontrol(...
    'Parent', uipanelId,...
    'Style', 'pushbutton',...
    'String', '',...
    'TooltipString', 'Apply multiplier',...
    'CData', applyIcon,...
    'BackgroundColor', gui.buttonBackgroundColor,...
    'Position', buttonPos(2, 1, gui.marginSub),...
    'CallBack',@applySel);  
  end


% __________________________ LAYERS FUNCTIONS _____________________________

  function drawLayerPanel(uipanelId)
    layUicontextmenu = uicontextmenu;
    uimenu(layUicontextmenu,...
      'Label', 'Line color',...
      'Callback', @changeSelLayColor);
    
    layStyleMenuId = uimenu(layUicontextmenu,...
      'Label', 'Line style');
    uimenu(layStyleMenuId,...
      'Label', 'Solid',...
      'Callback', {@changeSelLayStyle, '-'});
    uimenu(layStyleMenuId,...
      'Label', 'Dashed',...
      'Callback', {@changeSelLayStyle, '--'});
    uimenu(layStyleMenuId,...
      'Label', 'Dot',...
      'Callback', {@changeSelLayStyle, ':'});
    uimenu(layStyleMenuId,...
      'Label', 'Dash-dot',...
      'Callback', {@changeSelLayStyle, '-.'});
    uimenu(layStyleMenuId,...
      'Label', 'None',...
      'Callback', {@changeSelLayStyle, 'none'});
    
    layWidthMenuId = uimenu(layUicontextmenu,...
      'Label', 'Line width');
    uimenu(layWidthMenuId,...
      'Label', '1.0',...
      'Callback', {@changeSelLayWidth, 1});
    uimenu(layWidthMenuId,...
      'Label', '2.0',...
      'Callback', {@changeSelLayWidth, 2});
    uimenu(layWidthMenuId,...
      'Label', '3.0',...
      'Callback', {@changeSelLayWidth, 3});
    uimenu(layWidthMenuId,...
      'Label', '4.0',...
      'Callback', {@changeSelLayWidth, 4});
    uimenu(layWidthMenuId,...
      'Label', '5.0',...
      'Callback', {@changeSelLayWidth, 5});
    uimenu(layWidthMenuId,...
      'Label', '6.0',...
      'Callback', {@changeSelLayWidth, 6});
    uimenu(layWidthMenuId,...
      'Label', '7.0',...
      'Callback', {@changeSelLayWidth, 7});
    uimenu(layWidthMenuId,...
      'Label', '8.0',...
      'Callback', {@changeSelLayWidth, 8});
    uimenu(layWidthMenuId,...
      'Label', '9.0',...
      'Callback', {@changeSelLayWidth, 9});
    uimenu(layWidthMenuId,...
      'Label', '10.0',...
      'Callback', {@changeSelLayWidth, 10});

    layMarkerMenuId = uimenu(layUicontextmenu,...
      'Label', 'Marker');
    uimenu(layMarkerMenuId,...
      'Label', 'None',...
      'Callback', {@changeSelLayMarker, 'none'});    
    uimenu(layMarkerMenuId,...
      'Label', '+',...
      'Callback', {@changeSelLayMarker, '+'});
    uimenu(layMarkerMenuId,...
      'Label', 'o',...
      'Callback', {@changeSelLayMarker, 'o'});
    uimenu(layMarkerMenuId,...
      'Label', '*',...
      'Callback', {@changeSelLayMarker, '*'});
    uimenu(layMarkerMenuId,...
      'Label', '.',...
      'Callback', {@changeSelLayMarker, '.'});
    uimenu(layMarkerMenuId,...
      'Label', 'x',...
      'Callback', {@changeSelLayMarker, 'x'});
    uimenu(layMarkerMenuId,...
      'Label', 'Square',...
      'Callback', {@changeSelLayMarker, 'square'});
    uimenu(layMarkerMenuId,...
      'Label', 'Diamond',...
      'Callback', {@changeSelLayMarker, 'diamond'});
    uimenu(layMarkerMenuId,...
      'Label', '^',...
      'Callback', {@changeSelLayMarker, '^'});
    uimenu(layMarkerMenuId,...
      'Label', 'v',...
      'Callback', {@changeSelLayMarker, 'v'});
    uimenu(layMarkerMenuId,...
      'Label', '>',...
      'Callback', {@changeSelLayMarker, '>'});
    uimenu(layMarkerMenuId,...
      'Label', '<',...
      'Callback', {@changeSelLayMarker, '<'});
    uimenu(layMarkerMenuId,...
      'Label', 'Pentagram',...
      'Callback', {@changeSelLayMarker, 'pentagram'});
    uimenu(layMarkerMenuId,...
      'Label', 'Hexagram',...
      'Callback', {@changeSelLayMarker, 'hexagram'});
    
    
    layMarkerSizeMenuId = uimenu(layUicontextmenu,...
      'Label', 'Marker Size');
    uimenu(layMarkerSizeMenuId,...
      'Label', '2.0',...
      'Callback', {@changeSelLayMarkerSize, 2});
    uimenu(layMarkerSizeMenuId,...
      'Label', '3.0',...
      'Callback', {@changeSelLayMarkerSize, 3});
    uimenu(layMarkerSizeMenuId,...
      'Label', '4.0',...
      'Callback', {@changeSelLayMarkerSize, 4});
    uimenu(layMarkerSizeMenuId,...
      'Label', '5.0',...
      'Callback', {@changeSelLayMarkerSize, 5});
    uimenu(layMarkerSizeMenuId,...
      'Label', '6.0',...
      'Callback', {@changeSelLayMarkerSize, 6});
    uimenu(layMarkerSizeMenuId,...
      'Label', '7.0',...
      'Callback', {@changeSelLayMarkerSize, 7});
    uimenu(layMarkerSizeMenuId,...
      'Label', '8.0',...
      'Callback', {@changeSelLayMarkerSize, 8});
    uimenu(layMarkerSizeMenuId,...
      'Label', '9.0',...
      'Callback', {@changeSelLayMarkerSize, 9});
    uimenu(layMarkerSizeMenuId,...
      'Label', '10.0',...
      'Callback', {@changeSelLayMarkerSize, 10});
    
    uimenu(layUicontextmenu,...
      'Label', 'Edit label',...
      'Callback', @changeSelLayLabel);
    
    
    layTypeMenuId = uimenu(layUicontextmenu,...
      'Label', 'Layer type');
    uimenu(layTypeMenuId,...
      'Label', 'Constant gain',...
      'Callback', {@changeSelLayConvType, 'constGain'});
    uimenu(layTypeMenuId,...
      'Label', 'Gain with smoothed border',...
      'Callback', {@changeSelLayConvType, 'smoothBorder'});
    uimenu(layTypeMenuId,...
      'Label', 'Hole filling',...
      'Callback', {@changeSelLayConvType, 'fill'});
    uimenu(layTypeMenuId,...
      'Label', 'Hole filling with noise phase',...
      'Callback', {@changeSelLayConvType, 'fillNoise'});
    
    uimenu(layUicontextmenu,...
      'Label', 'Invert layer',...
      'Callback', @invertSelLay);
    
    uimenu(layUicontextmenu,...
      'Label', 'Duplicate layer',...
      'Callback', @duplicateSelLay);
    
    uimenu(layUicontextmenu,...
      'Label', 'Clear layer',...
      'Callback', @clearSelLay);
    
    uimenu(layUicontextmenu,...
      'Label', 'Delete layer',...
      'Callback', @delSelLay);
    
    uimenu(layUicontextmenu,...
      'Label', 'Add new layer',...
      'Callback', @addSelLay);
    
    gui.layListId = uicontrol(...
      'Parent', uipanelId,...
      'HandleVisibility', 'off',...
      'Style', 'listbox',...
      'Fontsize', gui.fontSize,...
      'String', '',...
      'UIContextMenu', layUicontextmenu,...
      'BackgroundColor', [1 1 1],...
      'Position', [gui.marginSub(1),...
        gui.marginSub(2)+2*gui.textHeight+2*gui.vertDist+1,...
        gui.panelWidth-2*gui.margin(1)-1, 60],...
      'CallBack',@changeSelLay);
    
    gui.layPanelId = uipanelId; 
    updateLayerList;
    changeSelLay;
  end

  function updateLayerList()
    layCell = {};
    if isfield(sel,'lay')
      for ind = 1:length(sel.lay)
        layCell{end+1} = sel.lay(ind).label;
      end
      set(gui.layListId, 'String', layCell);
    end
  end
  
  function changeSelLayColor(objId, eventData)
    sel.lay(sel.curLay).color = uisetcolor;
    delSel;
    drawAllSel;
  end

  function changeSelLayStyle(objId, eventData, lineStyle)
    sel.lay(sel.curLay).lineStyle = lineStyle;
    delSel;
    drawAllSel;
  end

  function changeSelLayWidth(objId, eventData, lineWidth)
    sel.lay(sel.curLay).lineWidth = lineWidth;
    delSel;
    drawAllSel;
  end

  function changeSelLayMarker(objId, eventData, marker)
    sel.lay(sel.curLay).marker = marker;
    delSel;
    drawAllSel;
  end

  function changeSelLayMarkerSize(objId, eventData, markerSize)
    sel.lay(sel.curLay).markerSize = markerSize;
    delSel;
    drawAllSel;
  end

  function changeSelLayLabel(objId, eventData)
    newLabel = inputdlg('New label for the current selection layer',...
      'Change layer label');
    sel.lay(sel.curLay).label = newLabel{1};
    updateLayerList;    
  end

  function changeSelLayConvType(objId, eventData, convType)
    delSel;
    sel.lay(sel.curLay).convType = convType;    
    
    % !!! see if I put the default value somewhere
    % and if I keep a memory of the preceeding values
    
    sel.lay(sel.curLay).param = [];
    
    switch convType
      case 'constGain'
        sel.lay(sel.curLay).param(1).name = 'Gain';
        sel.lay(sel.curLay).param(1).val = 0;
      case 'smoothBorder'
        sel.lay(sel.curLay).param(1).name = 'Gain';
        sel.lay(sel.curLay).param(1).val = 0;
        sel.lay(sel.curLay).param(2).name = 'Border';
        sel.lay(sel.curLay).param(2).val = 10;
      case 'fill'
        sel.lay(sel.curLay).param(1).name = 'Width';
        sel.lay(sel.curLay).param(1).val = 2;
        sel.lay(sel.curLay).param(2).name = 'Height';
        sel.lay(sel.curLay).param(2).val = 3;
      case 'fillNoise'
        sel.lay(sel.curLay).param(1).name = 'Width';
        sel.lay(sel.curLay).param(1).val = 2;
        sel.lay(sel.curLay).param(2).name = 'Height ';
        sel.lay(sel.curLay).param(2).val = 3;
    end
    changeSelLay;
    drawAllSel;    
  end

  function invertSelLay(objId, eventData)      
       
      if isempty(sel.lay(sel.curLay).poly)
        sel.lay(sel.curLay).poly = fullSigPoly;
      else
        % compute difference between a polygon around the whole signal and 
        % the current selection to inverse current selection
        newPoly = PolygonClip(fullSigPoly, sel.lay(sel.curLay).poly, 0);
        delCurPoly;
        sel.lay(sel.curLay).poly = newPoly;
      end
      
      drawCurSel;
  end

  function poly = fullSigPoly()
    % compute a polygon corresponding to the whole signal
    temp = convIndToAxesX([1, 2]);
    temp = (temp(2) - temp(1)) / 2;
    tempX = convIndToAxesX([0, size(coeff.oriC, 2)]);
    tempX = [tempX(1) - temp, tempX(2) + temp]; 
    temp = convIndToAxesY([1, 2]);
    temp = (temp(2) - temp(1)) / 2;
    tempY = convIndToAxesY([0, size(coeff.oriC, 1)]);
    tempY = [tempY(1) - temp, tempY(2) + temp]; 
    poly.x = [tempX(1); tempX(2); tempX(2); tempX(1)];
    poly.y = [tempY(1); tempY(1); tempY(2); tempY(2)];
    poly.hole = false;
  end

  function duplicateSelLay(objId, eventData)
    sel.lay(end+1) = sel.lay(sel.curLay);
    sel.lay(end).label = [sel.lay(end).label ' copy'];
    if isfield(sel.lay(end).poly, 'id')
      sel.lay(end).poly.id = [];
    end
    updateLayerList;
  end

  function clearSelLay(objId, eventData)    
    delSel;
    updateUndoData;
    sel.lay(sel.curLay).poly = [];
    drawAllSel;
  end

  function delSelLay(objId, eventData)
    % check that we leave at least one layer
    if length(sel.lay) == 1
      warndlg('The selection must contain at least one layer');
      return;
    end
    
    delSel;
    updateUndoData;
    sel.lay = sel.lay([1:sel.curLay-1 sel.curLay+1:end]);
    updateLayerList;
    set(gui.layListId, 'Value', 1);
    changeSelLay;
    drawAllSel;
  end

  function addSelLay(objId, eventData)
    sel.lay(end+1) = default.sel.lay;
    sel.lay(end).label = 'New layer';
    updateLayerList;
  end

  function changeSelLay(objId, eventData)
    sel.curLay = get(gui.layListId, 'Value');
    drawLayParam(sel.curLay);
  end

  function drawLayParam(layInd)
    if ~isfield(sel,'lay')
       return;
    end
    child = findobj(gui.layPanelId);
    child = setdiff(child, gui.layPanelId);
    delete(child);
    for ind = 1:length(sel.lay(layInd).param)
      uicontrol(...
        'Parent', gui.layPanelId,...
        'FontSize', gui.fontSize, ...
        'Style', 'text',...
        'String', sel.lay(layInd).param(ind).name,...
        'Position', [gui.marginSub(1),...
          gui.marginSub(2)+(ind-1)*(gui.textHeight+gui.vertDist),...
          gui.textWidth, gui.textHeight]) ;
      uicontrol(...
        'Parent', gui.layPanelId,...
        'FontSize', gui.fontSize, ...
        'Style', 'edit',...
        'String', num2str(sel.lay(layInd).param(ind).val),...
        'BackgroundColor', gui.buttonBackgroundColor,...
        'Position', [gui.marginSub(1)+gui.textWidth,...
          gui.marginSub(2)+(ind-1)*(gui.textHeight+gui.vertDist),...
          gui.editWidth, gui.textHeight],...
        'Callback', {@changeLayParam, layInd, ind});
    end
  end

  function changeLayParam(objId, eventData, layInd, paramInd)
    handleSymb;
    temp = str2num(get(objId, 'String'));
    if isempty(temp)
      errordlg([sel.lay(layInd).param(paramInd).name...
        ' parameter is invalid']);
      set(objId, 'String', num2str(sel.lay(layInd).param(paramInd).val));
      return;
    end
    sel.lay(layInd).param(paramInd).val = temp;
  end

  function changeOpacity(objId, eventData)
    temp = str2num(get(objId, 'String'));
    if isempty(temp) || temp < 0 || temp > 1
      errordlg(['Opacity parameter is invalid, it must be a number '...
        'between 0 (transparent) and 1 (opaque)']);
      set(objId, 'String', num2str(gui.symbOpacity));
      return;
    end
    gui.symbOpacity = temp;
    handleSymb();
  end

% __________________ VISUALIZATION TOOLS FUNCTIONS ________________________

  function drawVisualizationTool(uipanelId)
    subPanelHeight = gui.margin(2) + gui.buttonHeight +...
      gui.textHeight + gui.fontSize + 4;
   
    visPanelPos = get(uipanelId,'Position');
  
    symbolPanelHeight = 2*subPanelHeight;

     currPos = 0;
    symbolPanelId = uipanel(...
      'Parent', uipanelId,...
      'Title', 'Symbol',...
      'TitlePosition', 'centertop',...
      'FontSize', gui.fontSize, ...
      'Units', 'pixels',...
      'Position', [1, currPos+1, gui.panelWidth-4 , symbolPanelHeight]);
    currPos = currPos + symbolPanelHeight;
    
     gui.symbListId = uicontrol(...
      'Parent', uipanelId,...
      'HandleVisibility', 'off',...
      'Style', 'listbox',...
      'Fontsize', gui.fontSize,...
      'String', '',...
      'BackgroundColor', [1 1 1],...
      'Position', [gui.marginSub(1),...
        gui.marginSub(2)+2*gui.textHeight+2*gui.vertDist+6,...
        gui.panelWidth-2*gui.margin(1)-1, 50],...
      'CallBack',@changeSelSymb);
    

    % draw layer panel
    % drawLayerPanel(symbolPanelId);
   
    showSymbolIcon = imread([iconpath,'showsymbol.png']);
    symbId = uicontrol(...
      'Parent', symbolPanelId,...
      'Style', 'togglebutton',...
      'String', '',...
      'TooltipString', 'Show multiplier symbol',...
      'CData', showSymbolIcon,...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'Position', buttonPos(2, 1, gui.marginSub),...
      'CallBack',@clicSymb);

    uicontrol(...
      'Parent', symbolPanelId,...
      'FontSize', gui.fontSize, ...
      'Style', 'Text',...
      'FontSize', gui.fontSize, ...
      'String', 'Opacity',...
      'Position', [gui.marginSub(1),...
        gui.marginSub(2)+gui.buttonHeight,...
        gui.textWidth,...
        gui.textHeight],...
      'CallBack',@changeOpacity);
    
    opacityId = uicontrol(...
      'Parent', symbolPanelId,...
      'FontSize', gui.fontSize, ...
      'Style', 'edit',...
      'String', num2str(default.opacity),...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'Position', [gui.marginSub(1)+gui.textWidth,...
        gui.marginSub(2)+gui.buttonHeight,...
        gui.editWidth,...
        gui.textHeight],...
      'CallBack',@changeOpacity);
   
   
   
    
    colormapPanelId = uipanel(...
      'Parent', uipanelId,...
      'Title', 'Colormap',...
      'TitlePosition', 'centertop',...
      'FontSize', gui.fontSize, ...
      'Units', 'pixels',...
      'Position', [1, currPos+1,...
        gui.panelWidth-4 , subPanelHeight]);
     currPos = currPos + subPanelHeight;
    
    colormapIcon = imread([iconpath,'colormap.png']);
    
    buttonColormapId = uicontrol(...
      'Parent', colormapPanelId,...
      'Style', 'pushbutton',...
      'String', '',...
      'TooltipString', 'Edit colormap',...
      'CData', colormapIcon,...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'Position', buttonPos(2, 1, gui.marginSub),...
      'CallBack',@editColormap);
    
    uicontrol(...
      'Parent', colormapPanelId,...
      'FontSize', gui.fontSize, ...
      'Style', 'Text',...
      'FontSize', gui.fontSize, ...
      'String', 'Dynamic',...
      'Position', [gui.marginSub(1), gui.marginSub(2)+gui.buttonHeight,...
        gui.textWidth, gui.textHeight],...
      'CallBack',@changeOpacity);
    
    gui.editDynamicId = uicontrol(...
      'Parent', colormapPanelId,...
      'FontSize', gui.fontSize, ...
      'Style', 'edit',...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'String', num2str(default.dynamic),...
      'Position', [gui.marginSub(1)+gui.textWidth,...
        gui.marginSub(2)+gui.buttonHeight,...
        gui.editWidth,...
        gui.textHeight],...
      'CallBack',@changeDynamic);
    
    zoomPanelId = uipanel(...
      'Parent', uipanelId,...
      'Title', 'Zoom',...
      'TitlePosition', 'centertop',...
      'FontSize', gui.fontSize, ...
      'Units', 'pixels',...
      'Position', [1, currPos+1, subpanelSize(1)]);
    
    zoomInIcon = imread([iconpath,'zoomin.png']);
    
    buttonZoomInId = uicontrol(...
      'Parent', zoomPanelId,...
      'Style', 'togglebutton',...
      'String', '',...
      'TooltipString', 'Zoom in',...
      'CData', zoomInIcon,...
      'Min', false,...
      'Max', true,...
      'Value', false,...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'Position', buttonPos(1, 1, gui.marginSub),...
      'CallBack', {@changeTool, 'zoomIn'});
    
    gui.tool(end).buttonId = buttonZoomInId;
    gui.tool(end).name= 'zoomIn';
    gui.tool(end).function = @switchZoomIn;
    
    zoomOutIcon = imread([iconpath,'zoomout.png']);
    
    buttonZoomOutId = uicontrol(...
      'Parent', zoomPanelId,...
      'Style', 'togglebutton',...
      'String', '',...
      'TooltipString', 'Zoom out',...
      'CData', zoomOutIcon,...
      'Min', false,...
      'Max', true,...
      'Value', false,...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'Position', buttonPos(2, 1, gui.marginSub),...
      'CallBack', {@changeTool, 'zoomOut'});
    
    gui.tool(end+1).buttonId = buttonZoomOutId;
    gui.tool(end).name = 'zoomOut';
    gui.tool(end).function = @switchZoomOut;

    panIcon = imread([iconpath,'pan.png']);
    
    buttonPanId = uicontrol(...
      'Parent', zoomPanelId,...
      'Style', 'togglebutton',...
      'String', '',...
      'TooltipString', 'Pan',...
      'CData', panIcon,...
      'Min', false,...
      'Max', true,...
      'Value', false,...
      'BackgroundColor', gui.buttonBackgroundColor,...
      'Position', buttonPos(3, 1, gui.marginSub),...
      'CallBack', {@changeTool, 'pan'});
    
    gui.tool(end+1).buttonId = buttonPanId;
    gui.tool(end).name = 'pan';
    gui.tool(end).function = @switchPan;
    
  end

  function clicSymb(objId, eventData)
    button_state = get(objId,'Value');
    if button_state == get(objId,'Max')
       delSel;
       gui.showSymb = true;
       drawAllSel;
    elseif button_state == get(objId,'Min')
       delSel;
       gui.showSymb = false;
       drawAllSel;
    end 

  end

  function switchZoomIn(toolInd)
    % toolInd is not used but needed as I pass a parameter to this 
    % function !!! precise this
    if getMatlabVersion >= 7.3
      zoomId = zoom(gui.mainFigId);
      set(zoomId,...
        'Direction', 'in',...
        'Enable', 'on',...
        'ActionPostCallback', @updateExploreRect);
    else
      zoom(gui.mainFigId, 'on');
      % !!! this is undocumented in version 7.2 but it works, do something
      % more general
      zoom(gui.mainFigId, 'Direction', 'in'); 
      % !!! problem: we miss something here to force update rectangle in 
      % overview after a zoom
    end
  end

  function switchZoomOut(toolInd) 
    % toolInd is not used but needed as I pass a parameter to this function
    % !!! precise this
    if getMatlabVersion >= 7.3
      zoomId = zoom(gui.mainFigId);
      set(zoomId,...
        'Direction', 'out',...
        'Enable', 'on',...
        'ActionPostCallback', @updateExploreRect);
    else
      zoom(gui.mainFigId, 'on');
      % !!! this is undocumented in version 7.2 but it works, do something
      % more general
      zoom(gui.mainFigId, 'Direction', 'out');
      % !!! problem: we miss something here to force update rectangle in 
      % overview after a zoom
    end
  end
    
  function switchPan(toolInd)
    % toolInd is not used but needed as I pass a parameter to this function
    % !!! precise this

    if getMatlabVersion >= 7.3
      panId = pan(gui.mainFigId);
      set(panId,...
        'Enable', 'on',...
        'ActionPostCallback', @updateExploreRect);
    else
      pan(gui.mainFigId, 'on');
      % !!! problem: we miss something here to force update rectangle in 
      % overview after a pan
    end
  end

  function [res] = updateExploreRect(objId, eventData)
    oriInd = find(strcmp('originalMain', {visu.label}), 1);
    exploreSetRect(get(visu(oriInd).axesId, 'XLim'),...
      get(visu(oriInd).axesId, 'YLim'));
  end

  function editColormap(objId, eventData)
    colormapeditor(gui.mainFigId);
  end

  function changeDynamic(objId, eventData)
    oriInd = find(strcmp('originalMain', {visu.label}), 1);
    temp = str2num(get(gui.editDynamicId, 'String'));
    if isempty(temp) || temp < 0
      errordlg(['Dynamic parameter is invalid, it must be a positive '...
        'number in dB']);
      set(gui.editDynamicId, 'String', num2str(visu(oriInd).dynamic));
      return;
    end
    visu(oriInd).dynamic = temp;
    updateVisu(false, true);
  end

  function exploreOverview(axesId)
    oriInd = find(strcmp('originalMain', {visu.label}), 1);
    axLim = axis(visu(oriInd).axesId);
    
    xLim = axLim(1:2);
    yLim = axLim(3:4);
    
    indXData = [1, 2, 2, 1];
    indYData = [1, 1, 2, 2];
    
    gui.exploreRect.patchId = patch(...
      'Parent', axesId,...
      'XData', xLim(indXData),...
      'YData', yLim(indYData),...
      'FaceColor', 'none',...
      'LineWidth', 2,...
      'EdgeColor', 'w',... % !!! put the color and width as variables
      'ButtonDownFcn', @explorePatchButtonDown);
    
    gui.exploreRect.lineId = line(...
      'Parent', axesId,...
      'XData', xLim(indXData),...
      'YData', yLim(indYData),...
      'LineStyle', 'None',...
      'LineWidth', 2,...
      'Color', 'w',... % !!! put the color  and width as variables
      'Marker', 'o',...
      'MarkerSize', 7,...
      'ButtonDownFcn', @exploreLineButtonDown);
    
    patchId = gui.exploreRect.patchId;
    lineId = gui.exploreRect.lineId;
    figId = gui.mainFigId;
    initPoint = [];

    indX = [];
    indY = [];
    
    exploreAxesXLim = [];
    exploreAxesYLim = [];
    
    function exploreLineButtonDown(objId, eventData)
      % the coordinate that should be modified depends on the selected 
      % corner
      point = get(axesId,'CurrentPoint');
      
      xData = get(lineId, 'XData');
      yData = get(lineId, 'YData');
      xLim = xData(1:2);
      yLim = yData(2:3);
      
      [mi, indX] = min(abs(xLim - point(1,1)));
      [mi, indY] = min(abs(yLim - point(1,2)));
      
      set(figId,'WindowButtonMotionFcn',@exploreMotionLine);
      set(figId,'WindowButtonUpFcn',@exploreButtonUp);
    end
    
    function exploreMotionLine(objId, eventData)
      point = get(axesId,'CurrentPoint');
      xLim(indX) = point(1,1);
      yLim(indY) = point(1,2);
      
      % reorder the values
      xLimSort = sort(xLim);
      yLimSort = sort(yLim);
      
      set(patchId,...
        'XData', xLimSort(indXData),...
        'YData', yLimSort(indYData));
      set(lineId,...
        'XData', xLimSort(indXData),...
        'YData', yLimSort(indYData));
      oriInd = find(strcmp('originalMain', {visu.label}), 1);
      axis(visu(oriInd).axesId, [xLimSort, yLimSort]);
      % !!! this plot can be removed if it's to heavy
      drawnow;
    end
    
    
    function explorePatchButtonDown(objId, eventData)
      initPoint = get(axesId,'CurrentPoint');
      exploreAxesXLim = get(axesId, 'XLim');
      exploreAxesYLim = get(axesId, 'YLim');
      xData = get(patchId, 'XData');
      yData = get(patchId, 'YData');
      xLim = xData(1:2);
      yLim = yData(2:3);
      set(figId,'WindowButtonMotionFcn',@exploreMotionPatch);
      set(figId,'WindowButtonUpFcn',@exploreButtonUp);
    end
    
    
    function exploreMotionPatch(objId, eventData)
      point = get(axesId,'CurrentPoint');
      newXLim = xLim + (point(1,1)-initPoint(1,1));
      newYLim = yLim + (point(1,2)-initPoint(1,2));
      
      if ~(newXLim(2) < exploreAxesXLim(1) ||...
          newXLim(1) > exploreAxesXLim(2) ||...
          newYLim(2) < exploreAxesYLim(1) ||...
          newYLim(1) > exploreAxesYLim(2))
        set(patchId,...
          'XData', newXLim(indXData),...
          'YData', newYLim(indYData));
        set(lineId,...
          'XData', newXLim(indXData),...
          'YData', newYLim(indYData));
        oriInd = find(strcmp('originalMain', {visu.label}), 1);
        axis(visu(oriInd).axesId, [newXLim; newYLim]);
        % !!! this plot can be removed if it's to heavy
        drawnow;
      end
    end
    
    function exploreButtonUp(objId, eventData)
      set(figId,'WindowButtonMotionFcn','');
      set(figId,'WindowButtonUpFcn','');
      xData = get(patchId, 'XData');
      yData = get(patchId, 'YData');
      oriInd = find(strcmp('originalMain', {visu.label}), 1);
      axis(visu(oriInd).axesId, [xData(1:2); yData(2:3)]);      
    end
    
  end

  function exploreSetRect(xLim, yLim)
    indXData = [1, 2, 2, 1];
    indYData = [1, 1, 2, 2];
    if ishandle(gui.exploreRect.patchId)
      set(gui.exploreRect.patchId,...
        'XData', xLim(indXData),...
        'YData', yLim(indYData));
    end
    if ishandle(gui.exploreRect.lineId)
      set(gui.exploreRect.lineId,...
        'XData', xLim(indXData),...
        'YData', yLim(indYData));
    end
  end

  function changeSelSymb(objId, eventData)
    symbol.curSymb = get(gui.symbListId, 'Value');
    handleSymb();
  end

  function resetSymbol() % erase sel and put it to default value
    symbol = default.symbol;
    if isfield(gui, 'symbListId')
      if ishandle(gui.symbListId)
        set(gui.symbListId, 'Value', symbol.curSymb);
        drawSymbolParam(symbol.curSymb);
        updateSymbolList;
      end
    end
  end

  function updateSymbolList()
    symCell = {};
    if isfield(symbol,'data')
       for ind = 1:length(symbol.data)
         symCell{end+1} = symbol.data(ind).name;
       end
       set(gui.symbListId, 'String', symCell);
    end
  end

  function drawSymbolParam(symbInd)
  % nothing to do yet
  end

  function addSymbolListItem(symb)
    symbol.data(end+1).val = symb;
    symbol.data(end).name = 'Imported';
    updateSymbolList();
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





end

