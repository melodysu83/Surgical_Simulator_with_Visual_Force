function varargout = Surgical_Simulation(varargin)
global Start;
if ~exist('Start','var') || isempty(Start)
    Start = false;
end
% SURGICAL_SIMULATION MATLAB code for Surgical_Simulation.fig
%      SURGICAL_SIMULATION, by itself, creates a new SURGICAL_SIMULATION or raises the existing
%      singleton*.
%
%      H = SURGICAL_SIMULATION returns the handle to a new SURGICAL_SIMULATION or the handle to
%      the existing singleton*.
%
%      SURGICAL_SIMULATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SURGICAL_SIMULATION.M with the given input arguments.
%
%      SURGICAL_SIMULATION('Property','Value',...) creates a new SURGICAL_SIMULATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Surgical_Simulation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Surgical_Simulation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Surgical_Simulation

% Last Modified by GUIDE v2.5 05-Mar-2017 23:07:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Surgical_Simulation_OpeningFcn, ...
                   'gui_OutputFcn',  @Surgical_Simulation_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Functions I added manually.
function Plot_JoyStickBasisFnc(handles,run_level,ctrlx,ctrly)
    %axes(handles.axes4);
    persistent x_2 x_4 x_6 x_8 y_2 y_4 y_6 y_8;
    if isempty(x_2)
        ang = 0:0.01:2*pi; 
        x_2=0.25*cos(ang);
        y_2=0.25*sin(ang);
        x_4=0.5*cos(ang);
        y_4=0.5*sin(ang);
        x_6=0.75*cos(ang);
        y_6=0.75*sin(ang);
        x_8=cos(ang);
        y_8=sin(ang);
    end
    if run_level == 0
        plot(handles.axes4,x_2,y_2,'k',x_4,y_4,'k',x_6,y_6,'k',x_8,y_8,'k',0,0,'ko',[-1,1],[0,0],'k-',[0,0],[-1,1],'k-',[0,ctrlx],[0,ctrly],'b-',ctrlx,ctrly,'bo','MarkerFaceColor','b');
    else
        if run_level == 1
            plot(handles.axes4,x_2,y_2,'k',x_4,y_4,'k',x_6,y_6,'k',x_8,y_8,'k',[-1,1],[0,0],'k-',[0,0],[-1,1],'k-',ctrlx,ctrly,'bo','MarkerFaceColor','b');
        else
            plot(handles.axes4,x_2,y_2,'k',x_4,y_4,'k',x_6,y_6,'k',x_8,y_8,'k',[-1,1],[0,0],'k-',[0,0],[-1,1],'k-',ctrlx,ctrly,'ro','MarkerFaceColor','r');
        end
    end 

function Plot_ForceAnalysisFnc(handles,F)
    global Force_See ForceArray;
    global Force_Thres;
    if ~exist('Force_Thres','var') || isempty(Force_Thres)
        Force_Thres = 0.3;
    end  
    F_realscale = F*150;
    if Force_See
%         % Create 2 copies of y
%         bottomLine = F;
%         topLine = F;
%         % Set the values you don't want to get drawn to nan
%         bottomLine(F >= Force_Thres) = NaN; 
%         topLine(F <= Force_Thres) = NaN;
%         % Plot it out
%         plot(handles.axes3,0:0.01:10,F,'g-',0:0.01:10,bottomLine,'g-',0:0.01:10,topLine,'r-','LineWidth',2,'Markersize',5);
        plot(handles.axes3,[0 0],[0 0.4],'k-',0:0.01:10,F_realscale,'g-',[0 10],[Force_Thres Force_Thres],'r--');    
        set(handles.text31,'string',['Applied Force: ' num2str(round(F_realscale(1)*100000)/100000) ' (N)']);  
    else
        plot(handles.axes3,[0 0],[0 0.4],'k-');
        set(handles.text31,'string','Applied Force: Cannot See');
        set(handles.axes3,'Color',[0,0,0]);
    end
    ForceArray = round(F_realscale(1)*100000)/100000;

    

function Plot_SimulationFnc(hObject, eventdata, handles)
    % Todo: target, deformation max limit, force max limit(real value)
    global AreaArray CollisionArray DeformMaxAll;
    global terrain_poly1 terrain_poly2 ctrlx ctrly Start Refresh_tissue_and_tool Ko K1 K2 ko k1 k2 b1 b2 A B Gripper Force_See tissue_option;
    persistent tissue_terrainx tissue_terrainy1 tissue_terrainy2 tissue_deform1 tissue_deform2 tool_startx tool_starty tool_endx tool_endy difference1 difference2 RO_filter;
    persistent Xo_tissue1 X1_tissue1 X2_tissue1 Xo_tissue2 X1_tissue2 X2_tissue2 F_mag F_ori;
    persistent target_x target_y target_index Force_Thres;
    global tool_leverx tool_levery COMPLETE_STATUS JUST_MESSAGEBOX;
    persistent tool_endx1 tool_endy1 tool_endx2 tool_endy2 GRIPPER_LENGTH NEAR_TARGET_DISTANCE_THRESHOLD MY_SUCCESS_ICON MY_FAIL_ICON TOOL_STEP_SIZE;
    if isempty(tissue_terrainx) || Refresh_tissue_and_tool
        % initialize tissue
        RO_filter =  [(1:150),151,(150:-1:1)].^2/151^2;
        tissue_terrainx = linspace(0,10,1001); % total points 1001 (step size 0.01)
        x = 10*log(1:11)/log(11); % 11 points in the range of 0-10
        tissue_terrainy2=[];   
        switch tissue_option
            case 1
                % make sure upper tissue and lower tissue dont overlap
                while isempty(tissue_terrainy2) || ~isempty(find(tissue_terrainy2 - tissue_terrainy1 < 0.1, 1))...
                        || ~isempty(find(tissue_terrainy1 < 0, 1)) ||  ~isempty(find(tissue_terrainy2 > 10, 1))
                    tmp1 = 0.5*randi(3,1,11)+linspace(1,5,11);
                    tmp2 = 0.5*randi(3,1,11)+linspace(3.5,7.5,11);
                    y1 = min(tmp1,tmp2);
                    y2 = max(tmp1,tmp2);
                    terrain_poly1 = polyfit(x,y1,5);
                    terrain_poly2 = polyfit(x,y2,5); 
                    tissue_terrainy1 = polyval(terrain_poly1,tissue_terrainx);
                    tissue_terrainy2 = polyval(terrain_poly2,tissue_terrainx);
                end
                % initialize target
                target_index = [randi(400,1,1) randi(2,1,1)];
                target_x = tissue_terrainx(target_index(1));
             case 2
                 load('default_tissue_data.mat');
                 tissue_terrainy1 = default_tissue1_1;
                 tissue_terrainy2 = default_tissue1_2;
                 % initialize target
                 target_index = default_target1;
                 target_x = tissue_terrainx(target_index(1));
             case 3
                 load('default_tissue_data.mat');
                 tissue_terrainy1 = default_tissue2_1;
                 tissue_terrainy2 = default_tissue2_2;
                 % initialize target
                 target_index = default_target2;
                 target_x = tissue_terrainx(target_index(1));
             case 4
                 load('default_tissue_data.mat');
                 tissue_terrainy1 = default_tissue3_1;
                 tissue_terrainy2 = default_tissue3_2;
                 % initialize target
                 target_index = default_target3;
                 target_x = tissue_terrainx(target_index(1));
             case 5
                 load('default_tissue_data.mat');
                 tissue_terrainy1 = default_tissue4_1;
                 tissue_terrainy2 = default_tissue4_2;
                 % initialize target
                 target_index = default_target4;
                 target_x = tissue_terrainx(target_index(1));
             case 6
                 load('default_tissue_data.mat');
                 tissue_terrainy1 = default_tissue5_1;
                 tissue_terrainy2 = default_tissue5_2;
                 % initialize target
                 target_index = default_target5;
                 target_x = tissue_terrainx(target_index(1));
        end
        % initialize tissue deformation
        tissue_deform1 = tissue_terrainy1;
        tissue_deform2 = tissue_terrainy2;
        difference1 = zeros(1,1001);
        difference2 = zeros(1,1001);
        % initialize  tool
        tool_startx = 10;
        tool_starty = round((tissue_terrainy1(1001)+tissue_terrainy2(1001))*100/2.0)/100.0;
        tool_endx = tool_startx;
        tool_endy = tool_starty;
        Refresh_tissue_and_tool = false;
        Xo_tissue1 = 0;
        X1_tissue1 = 0;
        X2_tissue1 = 0; 
        Xo_tissue2 = 0; 
        X1_tissue2 = 0; 
        X2_tissue2 = 0;
        F_mag = zeros(1,1001);
        F_ori = 0;
        
        %initialize gripper
        Gripper = false;
        GRIPPER_LENGTH = 0.5;
        NEAR_TARGET_DISTANCE_THRESHOLD = 0.18;
        MY_SUCCESS_ICON = imread('success.png');
        MY_FAIL_ICON = imread('fail.png');
        TOOL_STEP_SIZE = 0.03;
        Force_Thres = 0.3;
        COMPLETE_STATUS = 0;
        % initialize statistics
        DeformMaxAll =0;
    end 
    
    % update tool (Todo: make more fancy)
        % (1) setup tool end point (virtual one, only true when gripper is closed)
        if ~isempty(Start)
            if Start && ~isempty(ctrlx) && ~isempty(ctrly)
                tool_endx = round(max(0,min(tool_endx + TOOL_STEP_SIZE*ctrlx,10))*100)/100.0; %round to 2 decimal
                tool_endy = round(max(0,min(tool_endy + TOOL_STEP_SIZE*ctrly,10))*100)/100.0;
                if tool_endx == 10 % avoid bug
                    tool_endx = tool_endx-0.01;
                end
            end
        end
        % (2) setup lever point
        if tool_startx - tool_endx > 0.01
            m = (tool_starty - tool_endy)/(tool_startx - tool_endx);
            k = tool_endy - m*tool_endx;
            deltax = GRIPPER_LENGTH/sqrt(1+m^2);
            tool_leverx = round(max(0,min(tool_endx + deltax,10))*100)/100.0;
            tool_levery = round(max(0,min(m*tool_leverx+k,10))*100)/100.0;
        else
             if tool_endy > tool_starty
                 tool_leverx = tool_startx;
                 tool_levery = round(max(tool_starty,tool_endy-GRIPPER_LENGTH)*100)/100.0;
             else
                 tool_leverx = tool_startx;
                 tool_levery = round(min(tool_endy+GRIPPER_LENGTH,tool_starty)*100)/100.0;
             end
        end
        % (3) setup tool end point (when gripper is open)
        if Gripper % Gripper open
            tmpx = (tool_endx - tool_leverx)*sqrt(3)/2+tool_leverx;
            tmpy = (tool_endy - tool_levery)*sqrt(3)/2+tool_levery;
            if exist('m','var') && ~isempty(m)
                theta1 = atan2(tool_endx-tool_startx,-tool_endy+tool_starty);   % for lower gripper half
                theta2 = atan2(-tool_endx+tool_startx,tool_endy-tool_starty);   % for upper gripper half   
                tool_endx1 =  round(max(0,min(tmpx - GRIPPER_LENGTH*cos(theta1)/2,10))*100)/100.0;    % lower gripper half
                tool_endy1 =  round(max(0,min(tmpy - GRIPPER_LENGTH*sin(theta1)/2,10))*100)/100.0;
                tool_endx2 =  round(max(0,min(tmpx - GRIPPER_LENGTH*cos(theta2)/2,10))*100)/100.0;    % upper gripper half
                tool_endy2 =  round(max(0,min(tmpy - GRIPPER_LENGTH*sin(theta2)/2,10))*100)/100.0;
            else
                tool_endx1 =  round(max(0,min(tmpx + GRIPPER_LENGTH/2,10))*100)/100.0;      % lower gripper half
                tool_endy1 =  round(max(0,min(tmpy,10))*100)/100.0;
                tool_endx2 =  round(max(0,min(tmpx - GRIPPER_LENGTH/2,10))*100)/100.0;      % upper gripper half
                tool_endy2 =  tool_endy1;
            end

        else % Gripper closed
            tool_endx1 = tool_endx;
            tool_endx2 = tool_endx;
            tool_endy1 = tool_endy;
            tool_endy2 = tool_endy;
        end
        
    % update tissue
        % (0) prepare for tissue update
        tooly1 = 10*ones(1,1001);
        tooly2 = -10*ones(1,1001);
        if tool_startx - tool_endx >0.01
            span = tool_leverx:0.01:tool_startx;
            tooly1(round(span*100+1)) = m * span + k;
            tooly2(round(span*100+1)) = m * span + k;
            if tool_leverx ~= tool_endx1
                m_1 = (tool_levery - tool_endy1)/(tool_leverx - tool_endx1);
                k_1 = tool_endy1 - m_1*tool_endx1;
                span1 = min(tool_endx1,tool_leverx):0.01:max(tool_endx1,tool_leverx);           
                tooly1(round(span1*100+1)) = min(m_1 * span1 + k_1,tooly1(round(span1*100+1)));
                tooly2(round(span1*100+1)) = max(m_1 * span1 + k_1,tooly2(round(span1*100+1)));
            else
                tooly1(round(tool_endx1*100+1)) = min(min(tool_endy1,tooly1(round(tool_endx1*100+1))),tool_levery);
                tooly2(round(tool_endx1*100+1)) = max(max(tool_endy1,tooly2(round(tool_endx1*100+1))),tool_levery);
            end
            if tool_leverx ~= tool_endx2
                m_2 = (tool_levery - tool_endy2)/(tool_leverx - tool_endx2);
                k_2 = tool_endy2 - m_2*tool_endx2;
                span2 = min(tool_endx2,tool_leverx):0.01:max(tool_endx2,tool_leverx);
                tooly1(round(span2*100+1)) = min(m_2 * span2 + k_2,tooly1(round(span2*100+1)));
                tooly2(round(span2*100+1)) = max(m_2 * span2 + k_2,tooly2(round(span2*100+1)));
            else
                tooly1(round(tool_endx2*100+1)) = min(min(tool_endy2,tooly1(round(tool_endx2*100+1))),tool_levery);
                tooly2(round(tool_endx2*100+1)) = max(max(tool_endy2,tooly2(round(tool_endx2*100+1))),tool_levery);
            end
            
        else
            if tool_endy > tool_starty % up right
                m_2 = (tool_levery - tool_endy2)/(tool_leverx - tool_endx2);
                k_2 = tool_endy2 - m_2*tool_endx2;
                span2 = min(tool_endx2,tool_leverx):0.01:max(tool_endx2,tool_leverx);
                tooly2(round(span2*100+1)) = m_2 * span2 + k_2;
            else % downward
                m_2 = (tool_levery - tool_endy2)/(tool_leverx - tool_endx2);
                k_2 = tool_endy2 - m_2*tool_endx2;
                span2 = min(tool_endx2,tool_leverx):0.01:max(tool_endx2,tool_leverx);
                tooly1(round(span2*100+1)) = m_2 * span2 + k_2;
            end
            tooly1(1001)=tool_endy;
            tooly2(1001)=tool_endy;
        end
        
        % (1) collision detection (Todo: make more fancy collision detection algorithm, and dividing space method)
        if tool_startx - tool_endx > 0.01
%           % method1 : all contact points
            check1_all = tissue_deform1-tooly1 >= 0;
            check2_all = tooly2-tissue_deform2 >= 0;
            % method2 : only the points start or ends intersection
            diff1 = diff(tissue_deform1-tooly1 >= 0);
            diff2 = diff(tooly2-tissue_deform2 >= 0);
            check1 = [0,diff1>0]+[diff1<0,0]~=0;
            check2 = [0,diff2>0]+[diff2<0,0]~=0;
            intersect_y = [tooly1(check1),tooly2(check2)];
            intersect_x = [tissue_terrainx(check1), tissue_terrainx(check2)];
            intersect_y_all = [tooly1(check1_all),tooly2(check2_all)];
            intersect_x_all = [tissue_terrainx(check1_all), tissue_terrainx(check2_all)];
        else
            if tool_endy < tissue_deform1(1001) || tool_endy > tissue_deform2(1001)
                intersect_y = tool_endy;
                intersect_x = tool_endx;
                check1 = zeros(1,1001);
                check2 = zeros(1,1001);
                check1(1001)=1;
                check2(1001)=1;
            end
        end
%         toolpose_x = tool_endx:tool_startx;
%         if tool_startx - tool_endx == 0 || length(toolpose_x) < 2
%             y1 = polyval(terrain_poly1,tool_endx);
%             y2 = polyval(terrain_poly2,tool_endx);
%             if tool_endy > y2
%                 intersect_x = tool_endx;
%                 intersect_y = y2;
%             else
%                 if tool_endy < y1
%                     intersect_x = tool_endx;
%                     intersect_y = y1;
%                 end
%             end        
%         else
%             m = (tool_starty - tool_endy)/(tool_startx - tool_endx);
%             k = tool_endy - m*tool_endx;
%             toolpose_y = m*toolpose_x+k;
%             [intersect_x1,intersect_y1] = intersections(toolpose_x,toolpose_y,tissue_terrainx,tissue_deform1,1);
%             [intersect_x2,intersect_y2] = intersections(toolpose_x,toolpose_y,tissue_terrainx,tissue_deform2,1);
%             intersect_x = [intersect_x1,intersect_x2];
%             intersect_y = [intersect_y1,intersect_y2];           
%         end      

        % (2) tissue deformation
        % stress relax from previous deform
        if exist('B','var') && ~isempty(B) && exist('difference1','var') && ~isempty(difference1) 
            small_def1 = difference1<0.01;
            small_def2 = difference2>-0.01;
            tissue_deform1(small_def1) = tissue_terrainy1(small_def1);
            tissue_deform2(small_def2) = tissue_terrainy2(small_def2);
            difference1(small_def1)=0;
            difference2(small_def2)=0;
            tissue_deform1 = tissue_deform1 + 200*B*sign(difference1).*abs(difference1).^A;
            tissue_deform2 = tissue_deform2 + 200*B*sign(difference2).*abs(difference2).^A;
        end
        % new deformation by intersection
        if exist('intersect_x','var')

            tissue_deform1 = min(tissue_deform1,tooly1);
            tissue_deform2 = max(tissue_deform2,tooly2);

            % Add edge smoothing for deformation
            difference1 = tissue_terrainy1-tissue_deform1; % >= 0
            difference2 = tissue_terrainy2-tissue_deform2; % <= 0
            intersect_index1 = find(check1 == 1);
            intersect_index2 = find(check2 == 1);
            for i = 1: length(intersect_index1)
                index = intersect_index1(i);
                tmp = zeros(1,1001);
                tmp(max(index-150,1):min(index+150,1001)) = difference1(index)*RO_filter(max(1,152-index):min(301,1152-index));
                difference1 = max(0,max(tmp,difference1));
            end
            for i = 1: length(intersect_index2)
                index = intersect_index2(i);
                tmp = zeros(1,1001);
                tmp(max(index-150,1):min(index+150,1001)) = difference2(index)*RO_filter(max(1,152-index):min(301,1152-index));
                difference2 = min(min(tmp,difference2),0);
            end
            tissue_deform1 =  min(tissue_terrainy1,min(10,max(0,tissue_terrainy1 - difference1)));
            tissue_deform2 =  max(tissue_terrainy2,min(10,max(0,tissue_terrainy2 - difference2)));
%             % method 2 : smoothing is not as good using this
%              tissue_deform1 = max(0,min(tissue_deform1,tissue_terrainy1-smooth(difference1,18,'lowess')'));
%              tissue_deform2 = min(10,max(tissue_deform2,tissue_terrainy2+smooth(-difference2,18,'lowess')'));            
        end
        
        % (3) force rendering  (Todo: horizontal force, quantify better, orientation)
        F_directionx = [];
        F_directiony = [];
        if exist('intersect_x','var') && ~isempty(intersect_x) 
            F_directionx = zeros(1,3*length(intersect_x));
            F_directiony = zeros(1,3*length(intersect_x));
            Area1 = trapz(difference1)/1000000.0;
            Area2 = trapz(-difference2)/1000000.0;
            AreaArray = max(Area1,Area2);
            if ~isempty(intersect_index1)  % intersect with lower tissue
                Xo_tissue1 = Area1;
                Xo1 = Xo_tissue1-X1_tissue1;
                Xo2 = Xo_tissue1-X2_tissue1;
                expXo1 = K1*Xo1*exp(k1*Xo1);
                expXo2 = K2*Xo2*exp(k2*Xo2);
                dX1 = min(Xo_tissue1-X1_tissue1,500*expXo1/b1);
                dX2 = min(Xo_tissue1-X2_tissue1,500*expXo2/b2);
                X1_tissue1 = X1_tissue1 + dX1;
                X2_tissue1 = X2_tissue1 + dX2;
                F_now = Ko*Xo_tissue1*exp(ko*Xo_tissue1)+expXo1+expXo2;
                if F_ori ~= 1
                    F_ori = 1;
                    Xo_tissue2 = 0;
                    X1_tissue2 = 0;
                    X2_tissue2 = 0;
                end
                theta = atan2(tool_leverx-tool_startx,-tool_levery+tool_starty);
                theta1 = atan2(tool_endx1-tool_leverx,-tool_endy1+tool_levery);
                theta2 = atan2(tool_endx2-tool_leverx,-tool_endy2+tool_levery);
                for i = 1:length(intersect_x) % Todo: orientation when there is only one tip point intersect
                    if intersect_x(i)==tool_endx1 && intersect_y(i)==tool_endy1
                        F_directionx(3*i-2:3*i) = [intersect_x(i) intersect_x(i)    NaN];
                        F_directiony(3*i-2:3*i)=  [intersect_y(i) intersect_y(i)+1  NaN];
                    else
                        if intersect_x(i)==tool_endx2 && intersect_y(i)==tool_endy2
                            F_directionx(3*i-2:3*i) = [intersect_x(i) intersect_x(i)    NaN];
                            F_directiony(3*i-2:3*i)=  [intersect_y(i) intersect_y(i)+1  NaN];
                        else
                            if exist('m','var') && ~isempty(m) && intersect_y(i) == m*intersect_x(i)+k
                                F_directionx(3*i-2:3*i)= [intersect_x(i) max(0,min(10,intersect_x(i)-cos(theta))) NaN];
                                F_directiony(3*i-2:3*i)= [intersect_y(i) max(0,min(10,intersect_y(i)-sin(theta))) NaN]; 
                            else
                                if Gripper && exist('m_1','var') && ~isempty(m_1) && intersect_y(i) == m_1*intersect_x(i)+k_1 
                                    F_directionx(3*i-2:3*i)= [intersect_x(i) max(0,min(10,intersect_x(i)-cos(theta1))) NaN];
                                    F_directiony(3*i-2:3*i)= [intersect_y(i) max(0,min(10,intersect_y(i)-sin(theta1))) NaN];

                                else 
                                    if Gripper && exist('m_2','var') && ~isempty(m_2) && intersect_y(i) == m_2*intersect_x(i)+k_2
                                        F_directionx(3*i-2:3*i)= [intersect_x(i) max(0,min(10,intersect_x(i)-cos(theta2))) NaN];
                                        F_directiony(3*i-2:3*i)= [intersect_y(i) max(0,min(10,intersect_y(i)-sin(theta2))) NaN];
                                    else
                                        F_directionx(3*i-2:3*i) = [intersect_x(i) max(0,min(10,intersect_x(i)))    NaN];
                                        F_directiony(3*i-2:3*i)=  [intersect_y(i) max(0,min(10,intersect_y(i)+1))  NaN];
                                    end
                                end
                            end
                        end
                    end 
                end
            else % intersect with upper tissue
                Xo_tissue2 = Area2;
                Xo1 = Xo_tissue2-X1_tissue2;
                Xo2 = Xo_tissue2-X2_tissue2;
                expXo1 = K1*Xo1*exp(k1*Xo1);
                expXo2 = K2*Xo2*exp(k2*Xo2);
                dX1 = min(Xo_tissue2-X1_tissue2,100*expXo1/b1);
                dX2 = min(Xo_tissue2-X2_tissue2,100*expXo2/b2);
                X1_tissue2 = X1_tissue2 + dX1;
                X2_tissue2 = X2_tissue2 + dX2;
                F_now = Ko*Xo_tissue2*exp(ko*Xo_tissue2)+expXo1+expXo2;
                if F_ori ~= 2
                    F_ori = 2;
                    Xo_tissue1 = 0;
                    X1_tissue1 = 0;
                    X2_tissue1 = 0;
                end
                theta = atan2(-tool_leverx+tool_startx,tool_levery-tool_starty);
                theta1 = atan2(-tool_endx1+tool_leverx,tool_endy1-tool_levery);
                theta2 = atan2(-tool_endx2+tool_leverx,tool_endy2-tool_levery);
                for i = 1:length(intersect_x) % Todo: orientation when there is only one tip point intersect
                    if intersect_x(i)==tool_endx1 && intersect_y(i)==tool_endy1
                        F_directionx(3*i-2:3*i) = [intersect_x(i) max(0,min(10,intersect_x(i)))    NaN];
                        F_directiony(3*i-2:3*i)=  [intersect_y(i) max(0,min(10,intersect_y(i)-1))  NaN];
                    else
                        if intersect_x(i)==tool_endx2 && intersect_y(i)==tool_endy2
                            F_directionx(3*i-2:3*i) = [intersect_x(i) max(0,min(10,intersect_x(i)))    NaN];
                            F_directiony(3*i-2:3*i)=  [intersect_y(i) max(0,min(10,intersect_y(i)-1))  NaN];
                        else
                            if exist('m','var') && ~isempty(m) && intersect_y(i) == m*intersect_x(i)+k
                                F_directionx(3*i-2:3*i)= [intersect_x(i) max(0,min(10,intersect_x(i)-cos(theta))) NaN];
                                F_directiony(3*i-2:3*i)= [intersect_y(i) max(0,min(10,intersect_y(i)-sin(theta))) NaN]; 
                            else
                                if Gripper && exist('m_1','var') && ~isempty(m_1) && intersect_y(i) == m_1*intersect_x(i)+k_1 
                                    F_directionx(3*i-2:3*i)= [intersect_x(i) max(0,min(10,intersect_x(i)-cos(theta1))) NaN];
                                    F_directiony(3*i-2:3*i)= [intersect_y(i) max(0,min(10,intersect_y(i)-sin(theta1))) NaN];

                                else 
                                    if Gripper && exist('m_2','var') && ~isempty(m_2) && intersect_y(i) == m_2*intersect_x(i)+k_2
                                        F_directionx(3*i-2:3*i)= [intersect_x(i) max(0,min(10,intersect_x(i)-cos(theta2))) NaN];
                                        F_directiony(3*i-2:3*i)= [intersect_y(i) max(0,min(10,intersect_y(i)-sin(theta2))) NaN];
                                    else
                                        F_directionx(3*i-2:3*i) = [intersect_x(i) max(0,min(10,intersect_x(i)))    NaN];
                                        F_directiony(3*i-2:3*i)=  [intersect_y(i) max(0,min(10,intersect_y(i)-1))  NaN];
                                    end
                                end
                            end
                        end
                    end 
                end
            end             
            F_mag = [max(0,F_now) F_mag(1:1000)];
            Plot_ForceAnalysisFnc(handles,F_mag);

%             dy = gradient(intersect_y);
%             dx = gradient(intersect_x);
%             normal_y = intersect_y+dx;
%             normal_x = intersect_x-dy;
        else
            F_now = 0;
            F_mag = [F_now F_mag(1:1000)];
            F_ori = 0;
            Plot_ForceAnalysisFnc(handles,F_mag);
        end

        % (4) optional: respiration (breathing : dynamic tissue)
        % Todo: add if we have time

    % update target
    if target_index(2) == 1  % on the lower 
        target_y = tissue_deform1(target_index(1))+0.05;
    else
        target_y = tissue_deform2(target_index(1))-0.05;
    end

    % Compute tool tip incision depth (Yana Code Merged)
%     intersect_index_at_tip = round(tool_endx*100)+1;
%     Distance1 = tool_endy - tissue_terrainy2(intersect_index_at_tip);
%     Distance2 = tissue_terrainy1(intersect_index_at_tip) - tool_endy;
%     if Distance1 > 0
%         set(handles.text33,'string',['Deform at Tip: ' num2str(Distance1) ' (cm)']); %output distance on the GUI
%     else
%         if Distance2 > 0
%             set(handles.text33,'string',['Deform at Tip: ' num2str(Distance2) ' (cm)']); %output distance on the GUI
%         else
%             set(handles.text33,'string',['Deform at Tip: ' num2str(0) ' (cm)']); %output distance on the GUI
%         end
%     end
    % Compute max incision depth
    deform_max1 = max(difference1);
    deform_max2 = max(-difference2);
    deform_max = max(deform_max1,deform_max2);
    if deform_max1 > deform_max2
        set(handles.text33,'string',['Max deform: ' num2str(deform_max) ' (cm)']); %output distance on the GUI
    else
        if deform_max1 < deform_max2
            set(handles.text33,'string',['Max deform: ' num2str(deform_max) ' (cm)']); %output distance on the GUI
        else
            set(handles.text33,'string',['Max deform: ' num2str(0) ' (cm)']); %output distance on the GUI
        end
    end    
    % prepare for statistics
    DeformMaxAll = max(DeformMaxAll,deform_max);
    
    % check if target is reached
    if tissue_option == 3 % tissue default 2 : the hardest , relax standard
        Distance_to_Target = (target_x - tool_endx)^2 + (target_y - tool_endy)^2;
        if Distance_to_Target < NEAR_TARGET_DISTANCE_THRESHOLD/10 && Gripper
            % prepare for statistics
            COMPLETE_STATUS = 1;
            Save_to_Statistics();
            % make sound
            [y,Fs] = audioread('gong.wav');
            sound(y,Fs); 
            % print message
            msgbox('Operation Completed!','Success','custom',MY_SUCCESS_ICON);
            pushbutton2_Callback(hObject, eventdata, handles);
            JUST_MESSAGEBOX = true;
        end
    else
        Distance_to_Target = (target_x - tool_leverx)^2 + (target_y - tool_levery)^2 + (target_x - tool_endx)^2 + (target_y - tool_endy)^2;
        In_Bound = tool_endx<target_x && target_x<tool_leverx && min(tool_endy,tool_levery)-0.01<=target_y && max(tool_endy,tool_levery)+0.01>=target_y;
        if Distance_to_Target < NEAR_TARGET_DISTANCE_THRESHOLD && Gripper && In_Bound
            % prepare for statistics
            COMPLETE_STATUS = 1;
            Save_to_Statistics();
            % make sound
            [y,Fs] = audioread('gong.wav');
            sound(y,Fs); 
            % print message
            msgbox('Operation Completed!','Success','custom',MY_SUCCESS_ICON);
            pushbutton2_Callback(hObject, eventdata, handles);
            JUST_MESSAGEBOX = true;
        end
    end   
    
    % check if exceed max force
    if F_mag(1)*150 > Force_Thres
        % prepare for statistics
        COMPLETE_STATUS = 0;
        Save_to_Statistics();
        % make sound
        [y,Fs] = audioread('gong.wav');
        sound(y,Fs); 
        % print message
        msgbox('Exceed Max Allowed Force!','Fail','custom',MY_FAIL_ICON);
        pushbutton2_Callback(hObject, eventdata, handles);
        JUST_MESSAGEBOX = true;
    end
     
    % plot tissue and tool
    if Gripper 	% Gripper open
        set(handles.text34,'string','Gripper:     Opened.');
        if ~exist('intersect_x','var')
            plot(handles.axes1,tissue_terrainx,tissue_terrainy1,'r:',tissue_terrainx,tissue_terrainy2,'r:',...
                tissue_terrainx,min(tissue_deform1,tissue_terrainy1),'r',tissue_terrainx,max(tissue_deform2,tissue_terrainy2),'r',[tool_startx, tool_leverx],...
                [tool_starty, tool_levery],'b-',[tool_endx1, tool_leverx],[tool_endy1, tool_levery],'b-',...
                [tool_endx2, tool_leverx],[tool_endy2, tool_levery],'b-',tool_startx,tool_starty,'.c',...
                [0,0],[0,10],'k',tool_leverx,tool_levery,'.c',target_x,target_y,'mo','MarkerFaceColor','r');
            plot(handles.axes2,[tool_startx, tool_leverx],[tool_starty, tool_levery],'b-',[tool_endx1, tool_leverx],...
                [tool_endy1, tool_levery],'b-',[tool_endx2, tool_leverx],[tool_endy2, tool_levery],'b-',tool_startx,...
                tool_starty,'.c',tool_leverx,tool_levery,'.c',[0,0],[0,10],'k',[0,10],[0,0],'k');
            set(handles.text8,'string','Collision Points: 0');
            CollisionArray = 0;
        else
            plot(handles.axes1,tissue_terrainx,tissue_terrainy1,'r:',tissue_terrainx,tissue_terrainy2,'r:',...
                tissue_terrainx,min(tissue_deform1,tissue_terrainy1),'r',tissue_terrainx,max(tissue_deform2,tissue_terrainy2),'r',[tool_startx, tool_leverx],...
                [tool_starty, tool_levery],'b-',[tool_endx1, tool_leverx],[tool_endy1, tool_levery],'b-',...
                [tool_endx2, tool_leverx],[tool_endy2, tool_levery],'b-',tool_startx,tool_starty,'.c',...
                [0,0],[0,10],'k',tool_leverx,tool_levery,'.c',target_x,target_y,'mo','MarkerFaceColor','r');
            if exist('F_directionx','var') && ~isempty(F_directionx)
                plot(handles.axes2,[tool_startx, tool_leverx],[tool_starty, tool_levery],'b-',[tool_endx1, tool_leverx],...
                    [tool_endy1, tool_levery],'b-',[tool_endx2, tool_leverx],[tool_endy2, tool_levery],'b-',tool_startx,...
                    tool_starty,'.c',tool_leverx,tool_levery,'.c',[0,0],[0,10],'k',[0,10],[0,0],'k',intersect_x_all,intersect_y_all,'g.',...
                    F_directionx, F_directiony,'r--');
            else
                plot(handles.axes2,[tool_startx, tool_leverx],[tool_starty, tool_levery],'b-',[tool_endx1, tool_leverx],...
                    [tool_endy1, tool_levery],'b-',[tool_endx2, tool_leverx],[tool_endy2, tool_levery],'b-',tool_startx,...
                    tool_starty,'.c',tool_leverx,tool_levery,'.c',[0,0],[0,10],'k',[0,10],[0,0],'k',intersect_x_all,intersect_y_all,'g');
            end   
            set(handles.text8,'string',['Collision Points: ' num2str(length(intersect_x_all))]);
            CollisionArray = length(intersect_x_all);
        end
    else        % Gripper closed
        set(handles.text34,'string','Gripper:     Closed.');
        if ~exist('intersect_x','var')
            plot(handles.axes1,tissue_terrainx,tissue_terrainy1,'r:',tissue_terrainx,tissue_terrainy2,'r:',...
                tissue_terrainx,tissue_deform1,'r',tissue_terrainx,tissue_deform2,'r',[tool_startx, tool_endx],...
                [tool_starty, tool_endy],'b-',tool_startx,tool_starty,'.c',...
                [0,0],[0,10],'k',tool_leverx,tool_levery,'.c',target_x,target_y,'mo','MarkerFaceColor','r');
            plot(handles.axes2,[tool_startx, tool_endx],[tool_starty, tool_endy],'b-',...
                                tool_startx,tool_starty,'.c',tool_leverx,tool_levery,'.c',[0,0],[0,10],'k',[0,10],[0,0],'k');
            set(handles.text8,'string','Collision Points: 0');
            CollisionArray = 0;
        else
            plot(handles.axes1,tissue_terrainx,tissue_terrainy1,'r:',tissue_terrainx,tissue_terrainy2,'r:',...
                tissue_terrainx,tissue_deform1,'r',tissue_terrainx,tissue_deform2,'r',[tool_startx, tool_endx],...
                [tool_starty, tool_endy],'b-',tool_startx,tool_starty,'.c',...
                [0,0],[0,10],'k',tool_leverx,tool_levery,'.c',target_x,target_y,'mo','MarkerFaceColor','r');
            if exist('F_directionx','var') && ~isempty(F_directionx)
                plot(handles.axes2,[tool_startx, tool_endx],[tool_starty, tool_endy],'b-',...
                                tool_startx,tool_starty,'.c',tool_leverx,tool_levery,'.c',[0,0],[0,10],'k',...
                                [0,10],[0,0],'k',intersect_x_all,intersect_y_all,'g.',F_directionx, F_directiony,'r--');
            else
                plot(handles.axes2,[tool_startx, tool_endx],[tool_starty, tool_endy],'b-',...
                                tool_startx,tool_starty,'.c',tool_leverx,tool_levery,'.c',...
                                [0,0],[0,10],'k',[0,10],[0,0],'k',intersect_x_all,intersect_y_all,'g');
            end   
            set(handles.text8,'string',['Collision Points: ' num2str(length(intersect_x_all))]);
            CollisionArray = length(intersect_x_all);
        end
    end
    
    % Optional: (Todo) color onlythe top and bottom region (not the middle)
    set(handles.axes1,'Color',[0.9255,0.8392,0.8392]);
    
    
    
function Set_TissueParamFcn(handles,command)
    global FirstPressStart;
    set(handles.popupmenu2,'enable',command,'Visible','on');
    set(handles.pushbutton1,'enable',command,'Visible','on');
    set(handles.slider1,'enable',command,'Visible','on');
    set(handles.slider2,'enable',command,'Visible','on');
    set(handles.slider3,'enable',command,'Visible','on');
    set(handles.slider4,'enable',command,'Visible','on');
    set(handles.slider5,'enable',command,'Visible','on');
    if FirstPressStart
        set(handles.pushbutton5,'enable',command,'Visible','on');
    else
        set(handles.pushbutton5,'enable','off','Visible','on');
    end

function Update_TissueParamFcn(handles)
    global K_slider k_slider b_slider A_slider B_slider Ko K1 K2 ko k1 k2 b1 b2 A B;
    global Ko_default K1_default K2_default ko_default k1_default k2_default b1_default b2_default A_default B_default;
    Ko = Ko_default*(1+0.1*K_slider);
    K1 = K1_default*(1+0.1*K_slider);
    K2 = K2_default*(1+0.1*K_slider);
    ko = ko_default*(1+0.1*k_slider);
    k1 = k1_default*(1+0.1*k_slider);
    k2 = k2_default*(1+0.1*k_slider);
    b1 = b1_default*(1+0.1*b_slider);
    b2 = b2_default*(1+0.1*b_slider);
    A = A_default*(1+0.5*A_slider);
    B = B_default*(1+0.1*B_slider);
    set(handles.slider1, 'value', K_slider);
    set(handles.slider2, 'value', k_slider);
    set(handles.slider3, 'value', b_slider);
    set(handles.slider4, 'value', A_slider);
    set(handles.slider5, 'value', B_slider);
    set(handles.text20,'string',['Ko = ' num2str(Ko)]);
    set(handles.text23,'string',['K1 = ' num2str(K1)]);
    set(handles.text24,'string',['K2 = ' num2str(K2)]);
    set(handles.text21,'string',['ko = ' num2str(ko)]);
    set(handles.text25,'string',['k1 = ' num2str(k1)]);
    set(handles.text26,'string',['k2 = ' num2str(k2)]);
    set(handles.text22,'string',['b1 = ' num2str(b1)]);
    set(handles.text27,'string',['b2 = ' num2str(b2)]);
    set(handles.text28,'string',['A = ' num2str(A)]);
    set(handles.text29,'string',['B = ' num2str(B)]);
    
function Update_TimerFnc(hObject, eventdata, handles)
global Start mytimer TOTAL_TIME PERIOD_TIME MY_FAIL_ICON ForceArray AreaArray CollisionArray tissue_option;
persistent filecount c file_name;
global ForceArrayAll CollisionArrayAll XArrayAll YArrayAll TimerArrayAll JUST_MESSAGEBOX;
global tool_leverx tool_levery COMPLETE_STATUS;
if ~exist('TOTAL_TIME','var') || isempty(TOTAL_TIME)
    TOTAL_TIME = 30;
    PERIOD_TIME = 0.033;
    MY_FAIL_ICON = imread('fail.png');
    filecount = 0;
    JUST_MESSAGEBOX = false;
end

if ~isempty(Start)
    if Start == true 
        if ~exist('filecount','var') || isempty(filecount) || filecount == -1 || JUST_MESSAGEBOX
            filecount = 0;
            c = clock;
            file_name = [num2str(c(1)) num2str(c(2)) num2str(c(3)) num2str(c(4)) num2str(c(5)) num2str(round(c(6))) '.txt'];
            
            TimerArrayAll = [];
            ForceArrayAll = [];
            CollisionArrayAll = [];
            XArrayAll = [];
            YArrayAll = [];
        end
        if ~exist('mytimer','var') || isempty(mytimer) || JUST_MESSAGEBOX
            mytimer = TOTAL_TIME;
            mysecond = round(mod(mytimer,60));
            mymicro = round((mytimer-round(mysecond))*100);
            if mymicro<0
                mymicro = 100+mymicro;
            end
            set(handles.text35,'string',['Timer Count Down:  ' num2str(mysecond) ' : ' num2str(mymicro)]);
            JUST_MESSAGEBOX = false;
        else
            mytimer = mytimer - PERIOD_TIME;           
            mysecond = round(mod(mytimer,60));
            mymicro = round((mytimer-mysecond)*100);
            if mymicro < 0
                mymicro = 100+mymicro;
            end
            if  mytimer < 0
                % complete status update
                COMPLETE_STATUS = 0;
                % show message
                msgbox('Operation timeout!','Fail','custom',MY_FAIL_ICON);
                pushbutton2_Callback(hObject, eventdata, handles);
                % prepare for statistics
                Save_to_Statistics();
                % make sound
                [y,Fs] = audioread('gong.wav');
                sound(y,Fs);
                mysecond = 0;
                mymicro = 0;
                JUST_MESSAGEBOX = true;
            end
            set(handles.text35,'string',['Timer Count Down:  ' num2str(mysecond) ' : ' num2str(mymicro)]);
            % write to file : for human subject experiment
            filecount = mod(filecount +1,15);
%             if filecount == 0
%                 FileArray = [mytimer ForceArray AreaArray CollisionArray tissue_option-1];
%                 dlmwrite(file_name,FileArray,'-append','delimiter',' ','roffset', 0,'newline','pc');
%             end

            % prepare for statistics
            TimerArrayAll = [TimerArrayAll TOTAL_TIME-mytimer];
            ForceArrayAll = [ForceArrayAll ForceArray];
            CollisionArrayAll = [CollisionArrayAll CollisionArray];
            XArrayAll = [XArrayAll tool_leverx];
            YArrayAll = [YArrayAll tool_levery];
        end
        Plot_SimulationFnc(hObject, eventdata, handles);
    else
        if c ~= -1
            mytimer = TOTAL_TIME;
            filecount = -1;
            c = -1;
        end
    end
end


% --- Executes just before Surgical_Simulation is made visible.
function Surgical_Simulation_OpeningFcn(hObject, eventdata, handles, varargin)
global PERIOD_TIME TOTAL_TIME MY_FAIL_ICON;
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Surgical_Simulation (see VARARGIN)

% Choose default command line output for Surgical_Simulation
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes Surgical_Simulation wait for user response (see UIRESUME)
% uiwait(handles.figure1);
if ~exist('TOTAL_TIME','var') || isempty(TOTAL_TIME)
    TOTAL_TIME = 30;
    PERIOD_TIME = 0.033;
    MY_FAIL_ICON = imread('fail.png');
end
time = timer('Period', PERIOD_TIME, 'ExecutionMode', 'fixedRate','TimerFcn', {@Update_TimerFnc, handles});
start(time);

% --- Outputs from this function are returned to the command line.
function varargout = Surgical_Simulation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
 

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global K_slider;
K_slider = get(hObject,'Value');
Update_TissueParamFcn(handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global k_slider;
k_slider = get(hObject,'Value');
Update_TissueParamFcn(handles);

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global b_slider;
b_slider = get(hObject,'Value');
Update_TissueParamFcn(handles);

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global A_slider;
A_slider = get(hObject,'Value');
Update_TissueParamFcn(handles);

% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global B_slider;
B_slider = get(hObject,'Value');
Update_TissueParamFcn(handles);

% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global K_slider k_slider b_slider A_slider B_slider;
K_slider = 0;
k_slider = 0;
b_slider = 0;
A_slider = 0;
B_slider = 0;
Update_TissueParamFcn(handles);



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global K_slider k_slider b_slider A_slider B_slider;
global Start Refresh_tissue_and_tool FirstPressStart JUST_MESSAGEBOX COMPLETE_STATUS;
FirstPressStart = true;
JUST_MESSAGEBOX = true;
Start = ~Start;
if isempty (Refresh_tissue_and_tool)
    Refresh_tissue_and_tool = false;
end
if Start
    str = 'Stop';  
    Set_TissueParamFcn(handles,'off');
    set(handles.text17,'string','Simulation starting. Push "Stop" to change tissue param.'); 
    Update_TissueParamFcn(handles);
    Refresh_tissue_and_tool = true;
else
    str = 'Start';
    Set_TissueParamFcn(handles,'on');
    set(handles.text17,'string','Push "Start" button when finish tuning tissue param.');
end
set(handles.pushbutton2,'string',str,'enable','on','Visible','on');

% --- Executes on mouse press over axes background.
function axes4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Start ctrlx ctrly;
persistent run_level;
if ~exist('Start','var') || isempty(Start)
    Start = false;
end
if Start
    pos = get(hObject,'currentpoint');
    x = pos(1); y = pos(2);
    
    % Mouse in joystick region:
    if y > 7 && y < 14 && x > 136 && x < 158
        run_level = 0;
        dy = 3.5;  % ymin 7, ymax 14
        dx = 11; % xmin 136, xmax 158
        midy = 10.5;
        midx = 147;
        ctrlx = min(max(-1,(x-midx)/dx),1);
        ctrly = min(max(-1,(y-midy)/dy),1);
        Plot_JoyStickBasisFnc(handles,run_level,ctrlx,ctrly);
    else
        if run_level ~= 1
            run_level = 1;
            ctrlx = 0;
            ctrly = 0;
            Plot_JoyStickBasisFnc(handles,run_level,ctrlx,ctrly);
        end
    end
else
    if run_level ~= 2
        run_level = 2;
        ctrlx = 0;
        ctrly = 0;
        Plot_JoyStickBasisFnc(handles,run_level,ctrlx,ctrly);
    end
end


% --- Executes during object creation, after setting all properties.
function pushbutton1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% print tissue model
imshow('tissue model1.png')
global K_slider k_slider b_slider A_slider B_slider;
global Ko_default K1_default K2_default ko_default k1_default k2_default b1_default b2_default A_default B_default;
Ko_default = 2.03;
K1_default = 0.438;
K2_default = 0.102;
ko_default = 909.9;
k1_default = 1522;
k2_default = 81.18;
b1_default = 5073;
b2_default = 39.24;
A_default = 1;
B_default = 0.0001;
K_slider = 0;
k_slider = 0;
b_slider = 0;
A_slider = 0;
B_slider = 0;


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
   global Gripper Force_See Start;
   % The keypressfcn for the figure.
   keyPressed = eventdata.Character;
    if exist('Start','var') && Start 
        if keyPressed == 'g' % Toggle Gripper open and close
            if ~exist('Gripper','var') || isempty(Gripper)
                Gripper = true;
            else
                Gripper = ~Gripper;
            end
        end
    else
        if exist('Start','var') && ~Start
            if keyPressed == 'f' % Toggle Force visualization
                if ~exist('Force_See','var') || isempty(Force_See)
                    Force_See = false;
                    plot(handles.axes3,[0 0],[0 0.4],'k-');
                    set(handles.text31,'string','Applied Force: Cannot See');
                    set(handles.axes3,'Color',[0,0,0]);
                else
                    Force_See = ~Force_See;
                    if Force_See
                        plot(handles.axes3,[0 0],[0 0.4],'w-');
                        set(handles.text31,'string','Applied Force: (N)');
                        set(handles.axes3,'Color',[1,1,1]);
                    else
                        plot(handles.axes3,[0 0],[0 0.4],'k-');
                        set(handles.text31,'string','Applied Force: Cannot See');
                        set(handles.axes3,'Color',[0,0,0]);
                    end
                end
            end
        end
    end


% --- Executes on key press with focus on pushbutton2 and none of its controls.
function pushbutton2_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
   global Gripper Force_See Start;
    % The keypressfcn for the figure.
    keyPressed = eventdata.Character;
    if exist('Start','var') && Start
        if keyPressed == 'g' % Toggle Gripper open and close
            if ~exist('Gripper','var') || isempty(Gripper)
                Gripper = true;
            else
                Gripper = ~Gripper;
            end
        end
    else
        if exist('Start','var') && ~Start
            if keyPressed == 'f' % Toggle Force visualization
                if ~exist('Force_See','var') || isempty(Force_See)
                    Force_See = false;
                    plot(handles.axes3,[0 0],[0 0.4],'k-');
                    set(handles.text31,'string','Applied Force: Cannot See');
                    set(handles.axes3,'Color',[0,0,0]);
                else
                    Force_See = ~Force_See;
                    if Force_See
                        plot(handles.axes3,[0 0],[0 0.4],'w-');
                        set(handles.text31,'string','Applied Force: (N)');
                        set(handles.axes3,'Color',[1,1,1]);
                    else
                        plot(handles.axes3,[0 0],[0 0.4],'k-');
                        set(handles.text31,'string','Applied Force: Cannot See');
                        set(handles.axes3,'Color',[0,0,0]);
                    end
                end
            end
        end
    end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
global tissue_option;
tissue_option = get(hObject,'Value');
        

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
global tissue_option;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
tissue_option = 1;

% --- Executes during object creation, after setting all properties.
function pushbutton5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global FirstPressStart;
FirstPressStart = false;

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% prepare for statistics
global tissue_option ForceArrayAll TimerArrayAll CollisionArrayAll XArrayAll YArrayAll DeformMaxAll COMPLETE_STATUS Force_See;
save('this_trial_data','TimerArrayAll','XArrayAll','YArrayAll',...
                        'ForceArrayAll','CollisionArrayAll','DeformMaxAll','tissue_option','Force_See','COMPLETE_STATUS');
% do the statistics
if exist('test_stats.m','file')
    run('test_stats.m');
else
    sprintf('[error]: missing statistics script.');
end

function Save_to_Statistics
global tissue_option ForceArrayAll TimerArrayAll DeformMaxAll COMPLETE_STATUS Force_See;
if exist('TimerArrayAll','var') && ~isempty(TimerArrayAll)
    if ~exist('stat_data.mat','file')
        tissue1_stat_data = [];
        tissue2_stat_data = [];
        tissue3_stat_data = [];
        tissue4_stat_data = [];
        tissue5_stat_data = [];
        tissuer_stat_data = [];
        tissue1_count = [0,0,0,0]; % force_success, force_all, no_force_success no_force_all
        tissue2_count = [0,0,0,0];
        tissue3_count = [0,0,0,0];
        tissue4_count = [0,0,0,0];
        tissue5_count = [0,0,0,0];
        tissuer_count = [0,0,0,0];
    else
        load('stat_data.mat');
    end
    switch tissue_option-1
            case 0 
                if Force_See
                    tissuer_count(1) = tissuer_count(1) + COMPLETE_STATUS;
                    tissuer_count(2) = tissuer_count(2) + 1;
                else
                    tissuer_count(3) = tissuer_count(3) + COMPLETE_STATUS;
                    tissuer_count(4) = tissuer_count(4) + 1;
                end
                tissuer_stat_data = [tissuer_stat_data;TimerArrayAll(length(TimerArrayAll)) max(ForceArrayAll) DeformMaxAll COMPLETE_STATUS];
            case 1
                if Force_See
                    tissue1_count(1) = tissue1_count(1) + COMPLETE_STATUS;
                    tissue1_count(2) = tissue1_count(2) + 1;
                else
                    tissue1_count(3) = tissue1_count(3) + COMPLETE_STATUS;
                    tissue1_count(4) = tissue1_count(4) + 1;
                end
                tissue1_stat_data = [tissue1_stat_data;TimerArrayAll(length(TimerArrayAll)) max(ForceArrayAll) DeformMaxAll COMPLETE_STATUS];
            case 2 
                if Force_See
                    tissue2_count(1) = tissue2_count(1) + COMPLETE_STATUS;
                    tissue2_count(2) = tissue2_count(2) + 1;
                else
                    tissue2_count(3) = tissue2_count(3) + COMPLETE_STATUS;
                    tissue2_count(4) = tissue2_count(4) + 1;
                end
                tissue2_stat_data = [tissue2_stat_data;TimerArrayAll(length(TimerArrayAll)) max(ForceArrayAll) DeformMaxAll COMPLETE_STATUS];
            case 3
                if Force_See
                    tissue3_count(1) = tissue3_count(1) + COMPLETE_STATUS;
                    tissue3_count(2) = tissue3_count(2) + 1;
                else
                    tissue3_count(3) = tissue3_count(3) + COMPLETE_STATUS;
                    tissue3_count(4) = tissue3_count(4) + 1;
                end
                tissue3_stat_data = [tissue3_stat_data;TimerArrayAll(length(TimerArrayAll)) max(ForceArrayAll) DeformMaxAll COMPLETE_STATUS];
            case 4 
                if Force_See
                    tissue4_count(1) = tissue4_count(1) + COMPLETE_STATUS;
                    tissue4_count(2) = tissue4_count(2) + 1;
                else
                    tissue4_count(3) = tissue4_count(3) + COMPLETE_STATUS;
                    tissue4_count(4) = tissue4_count(4) + 1;
                end
                tissue4_stat_data = [tissue4_stat_data;TimerArrayAll(length(TimerArrayAll)) max(ForceArrayAll) DeformMaxAll COMPLETE_STATUS];
            case 5
                if Force_See
                    tissue5_count(1) = tissue5_count(1) + COMPLETE_STATUS;
                    tissue5_count(2) = tissue5_count(2) + 1;
                else
                    tissue5_count(3) = tissue5_count(3) + COMPLETE_STATUS;
                    tissue5_count(4) = tissue5_count(4) + 1;
                end
                tissue5_stat_data = [tissue5_stat_data;TimerArrayAll(length(TimerArrayAll)) max(ForceArrayAll) DeformMaxAll COMPLETE_STATUS];
    end
     save('stat_data','tissue1_stat_data','tissue2_stat_data','tissue3_stat_data','tissue4_stat_data',...
                     'tissue5_stat_data','tissuer_stat_data','tissue5_count','tissue4_count','tissue3_count','tissue2_count','tissue1_count','tissuer_count');
end


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3
global Force_See;
Force_See = true;


% --- Executes on key press with focus on pushbutton5 and none of its controls.
function pushbutton5_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
   global Gripper Force_See Start;
    % The keypressfcn for the figure.
    keyPressed = eventdata.Character;
    if exist('Start','var') && Start
        if keyPressed == 'g' % Toggle Gripper open and close
            if ~exist('Gripper','var') || isempty(Gripper)
                Gripper = true;
            else
                Gripper = ~Gripper;
            end
        end
    else
        if exist('Start','var') && ~Start
            if keyPressed == 'f' % Toggle Force visualization
                if ~exist('Force_See','var') || isempty(Force_See)
                    Force_See = false;
                    plot(handles.axes3,[0 0],[0 0.4],'k-');
                    set(handles.text31,'string','Applied Force: Cannot See');
                    set(handles.axes3,'Color',[0,0,0]);
                else
                    Force_See = ~Force_See;
                    if Force_See
                        plot(handles.axes3,[0 0],[0 0.4],'w-');
                        set(handles.text31,'string','Applied Force: (N)');
                        set(handles.axes3,'Color',[1,1,1]);
                    else
                        plot(handles.axes3,[0 0],[0 0.4],'k-');
                        set(handles.text31,'string','Applied Force: Cannot See');
                        set(handles.axes3,'Color',[0,0,0]);
                    end
                end
            end
        end
    end
