function oh = CounterForInput (defaultVal,timeout)
defaultVal = 5;
timeout = 15;

f = figure('Name','Input','MenuBar','none','ToolBar','none',...
           'NumberTitle','off','Position',[500 500 300 120]);

txt = uicontrol('Style','text','Position',[50 70 200 40],...
                'String','Time left: 10 s','FontSize',12);

editBox = uicontrol('Style','edit','Position',[50 40 200 25],...
                    'String','');

btn = uicontrol('Style','pushbutton','Position',[100 10 100 25],...
                'String','OK','Callback','uiresume(gcbf)');

tStart = tic;

while toc(tStart) < timeout
    remaining = ceil(timeout - toc(tStart));
    set(txt, 'String', sprintf('Where to cut the contour?\n Time left: %d s', remaining));
    drawnow;
end

uiresume(f);  % continue if time runs out

userStr = get(editBox,'String');
close(f);

if isempty(userStr)
    oh = defaultVal;
    fprintf('No input received. Using default: %g\n', defaultVal);
else
    oh = str2double(userStr);
end
end