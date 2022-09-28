function [molecule,out]=ppfit_digest(px_to_kb, fr_to_sec)

%fr_to_sec is the same as sec/frame
% Made by Ilya Finkelstein - 12.16.13 for RecBCD
% Modified for use with Exo1 Resection by Logan Myler
%     * Updated to work with MATLAB R2013

[filenames, pathname] = uigetfile2('*.*','Select files to run through', 'MultiSelect','on');

numfiles=size(filenames,2);

for fnum=1:numfiles
    vals=importdata(fullfile(pathname,char(filenames(fnum))));  %load data from MAT or text file
    x=vals(1:end,1)';
    offset=vals(1:end,4)';
    y=vals(1:end,3)';
  
    %convert to real-world units
    x=x*fr_to_sec;
    y = sqrt(y.^2 + offset.^2);
    y=(mean(y(1:20))-y)*px_to_kb; %use first ~20 points as y-offset
    
    disp([num2str(fnum) ' of ' num2str(numfiles) ': ' char(filenames(fnum))]);
    
    %create molecule structure so that it can be populated later
    molecule(fnum).x=x;
    molecule(fnum).y=y;
    molecule(fnum).keep=true;   %keep this molecule in final analysis
    molecule(fnum).filename=char(filenames(fnum));
    molecule(fnum).number=fnum;
    molecule(fnum).start_length=0;
    molecule(fnum).end_length=0;
    molecule(fnum).digest_rate_2=[];
    figure(1);
    plot(x,y,'-b');
    title('a: done (keep trace); q: ignore file; s: fit start plateu; d: fit decay slope; x: second slope; f: fit final plateu');
    points_to_fit=[];
    k=0;
    while k==0
    %select points for exclusion
        k = waitforbuttonpress;
      
        if k==0 %mouse click
            point1 = get(gca,'CurrentPoint');    % button down detected
            finalRect = rbbox;                   % return figure units
            point2 = get(gca,'CurrentPoint');    % button up detected
            point1 = point1(1,1:2);              % extract x and y
            point2 = point2(1,1:2);

            %find all data points that fall within rectangle
            x_bla=find(x>min(point1(1), point2(1)) & x<max(point1(1), point2(1)));
            y_bla=find(y>min(point1(2), point2(2)) & y<max(point1(2), point2(2)));
            incl=intersect(x_bla, y_bla);
            points_to_fit=union(points_to_fit, incl);
            figure(1);
            plot(x,y,'-b',x(points_to_fit), y(points_to_fit), '-g')     %show selected points
            title('a: done (keep trace); q: ignore file; s: fit start plateu; d: fit decay slope; x: second slope; f: fit final plateu') 
        else %determine which button was pressed and do proper command 
            btn=get(gcf, 'CurrentCharacter');
            switch lower(btn)
                case 'q'    %drop this file, step out of while loop
                    molecule(fnum).keep=false;
                    break %while loop
                    
                case 's'   % fit selected data to a straight line to get initial DNA length
                    if length(points_to_fit)>2     %make sure there's some data already selected
                        pp=polyfit(x(points_to_fit), y(points_to_fit),0); %fit to a straight line
                        linefit=polyval(pp, x);
                        
                        %store fit results
                        molecule(fnum).start_length=pp;                    
                        molecule(fnum).start_points_to_fit=points_to_fit;
                        molecule(fnum).start_x=x(points_to_fit);
                        molecule(fnum).start_y=y(points_to_fit);
                        
                        
                        %show result
                        plot(x,y,'-b',x(points_to_fit), y(points_to_fit), '-g', x(points_to_fit), linefit(points_to_fit), '-r')     %show selected points
                        title('a: done (keep trace); q: ignore file; s: fit start plateu; d: fit decay slope; x: second slope; f: fit final plateu') 
                        points_to_fit=[];   %clear selected points
                    end %if
                case 'f'   % fit selected data to a straight line to get final DNA length
                    if length(points_to_fit)>2     %make sure there's some data already selected
                        pp=polyfit(x(points_to_fit), y(points_to_fit),0); %fit to a straight line
                        linefit=polyval(pp, x);
                        
                        %store fit results
                        molecule(fnum).end_length=pp;                    
                        molecule(fnum).end_points_to_fit=points_to_fit;
                        molecule(fnum).end_x=x(points_to_fit);
                        molecule(fnum).end_y=y(points_to_fit);
                        
                        %show result
                        plot(x,y,'-b',x(points_to_fit), y(points_to_fit), '-g', x(points_to_fit), linefit(points_to_fit), '-r')     %show selected points
                        title('a: done (keep trace); q: ignore file; s: fit start plateu; d: fit decay slope; x: second slope; f: fit final plateu') 
                        points_to_fit=[];   %clear selected points
                    end %if
                case 'd'   % fit selected data to a sloped line to get digestion speed
                    if length(points_to_fit)>2     %make sure there's some data already selected
                        pp=polyfit(x(points_to_fit), y(points_to_fit),1); %fit to a diagonal line
                        linefit=polyval(pp, x);
                        
                        %store fit results
                        molecule(fnum).digest_pp=pp;
                        molecule(fnum).digest_rate=-pp(1);
                        
                        molecule(fnum).digest_points_to_fit=points_to_fit;
                        molecule(fnum).digest_x=x(points_to_fit);
                        molecule(fnum).digest_y=y(points_to_fit);
                        
                        %show result
                        plot(x,y,'-b',x(points_to_fit), y(points_to_fit), '-g', x(points_to_fit), linefit(points_to_fit), '-r')     %show selected points
                        title('a: done (keep trace); q: ignore file; s: fit start plateu; d: fit decay slope; x: second slope; f: fit final plateu') 
                        points_to_fit=[];   %clear selected points
                    end %if
                case 'x'   % fit selected data to a sloped line to get different digestion speed
                    if length(points_to_fit)>2     %make sure there's some data already selected
                        pp=polyfit(x(points_to_fit), y(points_to_fit),1); %fit to a diagonal line
                        linefit=polyval(pp, x);
                        
                        %store fit results
                        molecule(fnum).digest_pp_2=pp;
                        molecule(fnum).digest_rate_2=-pp(1);
                        
                        molecule(fnum).digest_points_to_fit_2=points_to_fit;
                        molecule(fnum).digest_x_2=x(points_to_fit);
                        molecule(fnum).digest_y_2=y(points_to_fit);
                        
                        %show result
                        plot(x,y,'-b',x(points_to_fit), y(points_to_fit), '-c', x(points_to_fit), linefit(points_to_fit), '-r')     %show selected points
                        title('a: done (keep trace); q: ignore file; s: fit start plateu; d: fit decay slope; x: second slope; f: fit final plateu') 
                        points_to_fit=[];   %clear selected points
                    end %if                    
                case 'a'    %done; keep trace
                    if molecule(fnum).keep==true
                        [path,name,ext]=fileparts(char(filenames(fnum)));
                        save_filename=   fullfile(pathname,[name '_analyzed.mat']);
                        molec=molecule(fnum);
                        if molecule(fnum).digest_rate_2==[]
                            molecule(fnum).digest_rate_2=0;
                        end
                        save(save_filename, 'molec');
                    end
                    break
                otherwise
                title('a: done (keep trace); q: ignore file; s: fit start plateu; d: fit decay slope; x: after chi slope; f: fit final plateu')            
            end %switch           
            k=0;
        end %k==0 if

    end %while loop
    
end %for loop over all files

counter=0;
out=[];
for fnum=1:numfiles        %loop through all files again to prepare output matrix
    if molecule(fnum).keep==true
        counter=counter+1;
            out(counter, 1)=molecule(fnum).number;
            out(counter, 2)=molecule(fnum).start_length;
            out(counter, 3)=molecule(fnum).end_length;        
            out(counter, 4)=-(molecule(fnum).start_length-molecule(fnum).end_length);   
            out(counter, 5)=-molecule(fnum).digest_rate;  
            out(counter, 6)=-molecule(fnum).digest_rate_2
    end
end
    
