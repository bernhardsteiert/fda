function ex_uicontrol
    % Example code for uicontrol reference page
    % http://www.mathworks.de/de/help/matlab/ref/uicontrol.html

    % Create a figure and an axes to 
    % contain a 3-D surface plot.
    figure
    hax = axes('Units','pixels');
    myextension = '130722_corrected';
    global angles
    angles = [0 0 0];

    load(sprintf('harm_basis_%s',myextension))
    time_range = getbasisrange(harm_basis);
    time_eval = linspace(time_range(1),time_range(2),201);
    basis_eval = eval_basis(harm_basis,time_eval);
    plot(time_eval,basis_eval,time_eval,time_eval*0,'k:')
    text(150,.18,sprintf('Angle 1 = %g',angles(1)))
    text(150,.16,sprintf('Angle 2 = %g',angles(2)))
    text(150,.14,sprintf('Angle 3 = %g',angles(3)))
    set(gca,'YLim',[-.2 .25])

    % Add a slider uicontrol to control the vertical 
    % scaling of the surface object. Position it 
    % under the Clear button.
    hSlider1 = uicontrol('Style', 'slider',...
        'Min',0,'Max',360,'Value',0,...
        'Position', [100 20 120 20],...
        'Callback', {@plot_harm_rot1,time_eval,basis_eval}); 
					% Slider function handle callback
          % Implemented as a local function
    hListener1 = handle.listener(hSlider1,'ActionEvent',{@plot_harm_rot1,time_eval,basis_eval});
    setappdata(hSlider1,'sliderListener',hListener1);
    hSlider2 = uicontrol('Style', 'slider',...
        'Min',0,'Max',360,'Value',0,...
        'Position', [250 20 120 20],...
        'Callback', {@plot_harm_rot2,time_eval,basis_eval}); 
					% Slider function handle callback
          % Implemented as a local function
    hListener2 = handle.listener(hSlider2,'ActionEvent',{@plot_harm_rot2,time_eval,basis_eval});
    setappdata(hSlider2,'sliderListener',hListener2);
    hSlider3 = uicontrol('Style', 'slider',...
        'Min',0,'Max',360,'Value',0,...
        'Position', [400 20 120 20],...
        'Callback', {@plot_harm_rot3,time_eval,basis_eval}); 
					% Slider function handle callback
          % Implemented as a local function
    hListener3 = handle.listener(hSlider3,'ActionEvent',{@plot_harm_rot3,time_eval,basis_eval});
    setappdata(hSlider3,'sliderListener',hListener3);
   
    % Add a text uicontrol to label the slider.
    uicontrol('Style','text',...
        'Position',[400 45 120 20],...
        'String','Angle 3')
    uicontrol('Style','text',...
        'Position',[250 45 120 20],...
        'String','Angle 2')
    uicontrol('Style','text',...
        'Position',[100 45 120 20],...
        'String','Angle 1')
end


function plot_harm_rot1(hObj,event,time_eval,basis_eval) %#ok<INUSL>
    % Called to set zlim of surface in figure axes
    % when user moves the slider control 
%     val = 51 - get(hObj,'Value');
    global angles
    
    basis_eval_rot = nan(size(basis_eval));
    angle1 = get(hObj,'Value');
    angles(1) = angle1;
    angle2 = angles(2);
    angle3 = angles(3);
    Rmat1 = [cos(2*pi*angle1/360)   sin(2*pi*angle1/360)   0; ...
             -sin(2*pi*angle1/360)  cos(2*pi*angle1/360)   0; ...
             0                      0                      1];
    Rmat2 = [cos(2*pi*angle2/360)   0  sin(2*pi*angle2/360); ...
             0                      1  0; ...
             -sin(2*pi*angle2/360)  0  cos(2*pi*angle2/360)];
    Rmat3 = [1 0                      0; ...
             0 cos(2*pi*angle3/360)   sin(2*pi*angle3/360); ...
             0 -sin(2*pi*angle3/360)  cos(2*pi*angle3/360)];
    Rmat = Rmat3 * Rmat2 * Rmat1;
    for irot = 1:size(Rmat,1)
        basis_eval_rot(:,irot) = sum(repmat(Rmat(irot,:),size(basis_eval,1),1).*basis_eval,2);
    end
    plot(time_eval,basis_eval_rot,time_eval,time_eval*0,'k:')
    text(150,.18,sprintf('Angle 1 = %g',angles(1)))
    text(150,.16,sprintf('Angle 2 = %g',angles(2)))
    text(150,.14,sprintf('Angle 3 = %g',angles(3)))
    set(gca,'YLim',[-.2 .25])
end

function plot_harm_rot2(hObj,event,time_eval,basis_eval) %#ok<INUSL>
    % Called to set zlim of surface in figure axes
    % when user moves the slider control 
%     val = 51 - get(hObj,'Value');
    global angles
    
    basis_eval_rot = nan(size(basis_eval));
    angle2 = get(hObj,'Value');
    angles(2) = angle2;
    angle1 = angles(1);
    angle3 = angles(3);
    Rmat1 = [cos(2*pi*angle1/360)   sin(2*pi*angle1/360)   0; ...
             -sin(2*pi*angle1/360)  cos(2*pi*angle1/360)   0; ...
             0                      0                      1];
    Rmat2 = [cos(2*pi*angle2/360)   0  sin(2*pi*angle2/360); ...
             0                      1  0; ...
             -sin(2*pi*angle2/360)  0  cos(2*pi*angle2/360)];
    Rmat3 = [1 0                      0; ...
             0 cos(2*pi*angle3/360)   sin(2*pi*angle3/360); ...
             0 -sin(2*pi*angle3/360)  cos(2*pi*angle3/360)];
    Rmat = Rmat3 * Rmat2 * Rmat1;
    for irot = 1:size(Rmat,1)
        basis_eval_rot(:,irot) = sum(repmat(Rmat(irot,:),size(basis_eval,1),1).*basis_eval,2);
    end
    plot(time_eval,basis_eval_rot,time_eval,time_eval*0,'k:')
    text(150,.18,sprintf('Angle 1 = %g',angles(1)))
    text(150,.16,sprintf('Angle 2 = %g',angles(2)))
    text(150,.14,sprintf('Angle 3 = %g',angles(3)))
    set(gca,'YLim',[-.2 .25])
end

function plot_harm_rot3(hObj,event,time_eval,basis_eval) %#ok<INUSL>
    % Called to set zlim of surface in figure axes
    % when user moves the slider control 
%     val = 51 - get(hObj,'Value');
    global angles
    
    basis_eval_rot = nan(size(basis_eval));
    angle3 = get(hObj,'Value');
    angles(3) = angle3;
    angle2 = angles(2);
    angle1 = angles(1);
    Rmat1 = [cos(2*pi*angle1/360)   sin(2*pi*angle1/360)   0; ...
             -sin(2*pi*angle1/360)  cos(2*pi*angle1/360)   0; ...
             0                      0                      1];
    Rmat2 = [cos(2*pi*angle2/360)   0  sin(2*pi*angle2/360); ...
             0                      1  0; ...
             -sin(2*pi*angle2/360)  0  cos(2*pi*angle2/360)];
    Rmat3 = [1 0                      0; ...
             0 cos(2*pi*angle3/360)   sin(2*pi*angle3/360); ...
             0 -sin(2*pi*angle3/360)  cos(2*pi*angle3/360)];
    Rmat = Rmat3 * Rmat2 * Rmat1;
    for irot = 1:size(Rmat,1)
        basis_eval_rot(:,irot) = sum(repmat(Rmat(irot,:),size(basis_eval,1),1).*basis_eval,2);
    end
    plot(time_eval,basis_eval_rot,time_eval,time_eval*0,'k:')
    text(150,.18,sprintf('Angle 1 = %g',angles(1)))
    text(150,.16,sprintf('Angle 2 = %g',angles(2)))
    text(150,.14,sprintf('Angle 3 = %g',angles(3)))
    set(gca,'YLim',[-.2 .25])
end
