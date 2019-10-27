
function[]=Project2
close all
clc
%Variables
%Figure variables
WALL_WIDTH = 3;
WALL_COLOR = [.3, .3, .8]; %format string for drawing walls
FIGURE_COLOR = [0, 0, 0]; %program background
AXIS_COLOR = [.15, .15, .15]; %the court
FIGURE_WIDTH = 900; %pixels
FIGURE_HEIGHT = 580;
PLOT_W = 150; %width in plot units. this will be main units for program
PLOT_H = 100; %height
GOAL_SIZE = 50;
GOAL_TOP = (PLOT_H+GOAL_SIZE)/2;
GOAL_BOT = (PLOT_H-GOAL_SIZE)/2;
CENTERLINE_COLOR =  [.2, .2, .8]; %format string for centerline
%ball variables
ballPlot = []; %main plot, includes ball and walls
ballX = []; %ball location
ballY = [];
BALL_MARKER_SIZE = 10; %aesthetic, does not affect physics, see BALL_RADIUS
BALL_COLOR = [.1, .7, .1];
BALL_OUTLINE = [.7, 1, .7];
BALL_SHAPE = 'o';
MIN_BALL_SPEED = .8;
FRAME_DELAY= 0.01;
BALL_RADIUS = 1.5; %radius to calculate bouncing
coll = 0;
X2 = 1:1:500;
Y2 = 3*X2+3;
PolyX = 20:1:40; %polygone X coordinate
PolyY = 0:1:40; % polygone Y coordinate


 function createFigure

%Drawing the screen and path
    %ScreenSize is a four-element vector: [left, bottom, width, height]:
    scrsz = get(0,'ScreenSize');
    fig = figure('Position',[(scrsz(3)-FIGURE_WIDTH)/2 ...
      (scrsz(4)-FIGURE_HEIGHT)/2 ...
      FIGURE_WIDTH, FIGURE_HEIGHT]);
    %register keydown and keyup listeners
    set(fig,'KeyPressFcn',@keyDown);
    %figure can't be resized
    set(fig, 'Resize', 'off');
    axis([0 PLOT_W 0 PLOT_H]);
    axis manual;
    %set color for the court, hide axis ticks.
    set(gca, 'color', AXIS_COLOR, 'YTick', [], 'XTick', []);
    %set background color for figure
    set(fig, 'color', FIGURE_COLOR);
    hold on;
    %plot walls
    topWallXs = [0,0,PLOT_W,PLOT_W];
    topWallYs = [GOAL_TOP,PLOT_H,PLOT_H,GOAL_TOP];
    bottomWallXs = [0,0,PLOT_W,PLOT_W];
    bottomWallYs = [GOAL_BOT,0,0,GOAL_BOT];
    plot(topWallXs, topWallYs, '-', ...
      'LineWidth', WALL_WIDTH, 'Color', WALL_COLOR);
    plot(bottomWallXs, bottomWallYs, '-', ...
      'LineWidth', WALL_WIDTH, 'Color', WALL_COLOR);
       centerline = plot([PLOT_W/2, PLOT_W/2],[PLOT_H, 0],'--');
        set(centerline, 'Color', CENTERLINE_COLOR);
        %Lines
        X = [10 200 799];
        Y = [1 300 479];
        Linia1 = line(X,Y,'Color','y','LineWidth',4)
         X1 = [100 500];
         Y1 = [1 799];
        Linia2 = line(X1,Y1,'Color','r','LineWidth',4)
        X2 = 1:1:500;
        Y2 = 3*X2+3;
        Linia3 = line (X2,Y2, 'Color', 'm', 'LineWidth', 1)
         ballPlot = plot(10,10);
         %Polygone
         X3 = [20 40 40 20];
         Y3 = [0 0 20 20];
         Polygone = patch (X3, Y3, [1 0 0], 'EdgeColor','green','FaceColor','none','LineWidth',2)
         %Sine wave 
         SineWave =10*sin(X2)+40; %Sine wave with offset
         plot (SineWave,X2)
    set(ballPlot, 'Marker', BALL_SHAPE);
    set(ballPlot, 'MarkerEdgeColor', BALL_OUTLINE);
    set(ballPlot, 'MarkerFaceColor', BALL_COLOR);
    set(ballPlot, 'MarkerSize', BALL_MARKER_SIZE);
 end

function moveBall
   ballX=0;
   ballY=0;
%  newX = ballX + (ballSpeed * ballVector(1));
%  newY = ballY + (ballSpeed * ballVector(2));
%  %Plot the ball
%  ballX = newX;
%  ballY = newY;
end

function refreshPlot
    set(ballPlot, 'XData', ballX, 'YData', ballY);
    drawnow;
    pause(FRAME_DELAY);
end
%------------keyDown------------
%listener registered in createFigure
%listens for input
%sets appropriate variables and calls functions
  function keyDown(src,event)
    switch event.Key
      case 'uparrow'
        ballX = ballX;
        ballY = ballY+1;
        refreshPlot;

      case 'downarrow'
        ballX = ballX;
        ballY = ballY-1;
        refreshPlot;
        
       case 'leftarrow'
          ballX = ballX-1;
          ballY = ballY;
          refreshPlot;
          
       case 'rightarrow'
          ballX = ballX+1;
          ballY = ballY;
          refreshPlot;
    end
    coll = collision (ballX, ballY, coll, X2, Y2, PolyX, PolyY);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Collision Detection Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [coll] = collision (ballX, ballY, coll, X2, Y2, PolyX, PolyY)
  L = length (X2);
  PolyX_max = max(PolyX);
  PolyX_min = min(PolyX);
  PolyY_max = max(PolyY);
  PolyY_min = min(PolyY);
  
  for k = 1:L
     if ballX == X2(k) && ballY == Y2(k)
        coll = coll+1;
     else
        coll = coll;
     end

  end
 if ballX > PolyX_min && ballX < PolyX_max && ballY > PolyY_min && ballY < PolyY_max
     coll = coll+1;
 else
     coll=coll;
 end
     disp (coll)
  end

createFigure;
moveBall;
refreshPlot;


end