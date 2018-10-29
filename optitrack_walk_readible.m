%% Taking data from Arena
%
% BVH is a text file which contains skeletal data, but its contents needs
% additional processing to draw the wireframe and create the animation.

name = '02_01_01';
[original_skeleton,time2] = loadbvh(name);
[compare_with_not_opt_skel,time4] = loadbvh(name);
[compare_with_opt_skel,time] = loadbvh(name);
[optimized_skel,time1] = loadbvh_2(name);
[not_optimized_skel,time3] = loadbvh_2(name);

%%
 global all_pt;
    
    global p1;
    global p2;
    global p3;
write_video = false;
all_angle=[0 0;];

% Prepare the new video file.
if write_video, vidObj = VideoWriter(name); open(vidObj); end
   
fincr = 10; 
for ff = 1:fincr:length(time) 
    h = figure(1); clf; hold on
    title(sprintf('%1.2f seconds',time(ff)))
    set(h,'color','white')
    
    %%
    %MAKE SKELETONS ON TO A SAME LEVEL
    for nn = 1:length(compare_with_opt_skel) 
        optimized_skel(nn).Dxyz(2,ff)=optimized_skel(nn).Dxyz(2,ff)+15; 
    end
    for nn = 1:length(compare_with_opt_skel) 
        original_skeleton(nn).Dxyz(1,ff)=original_skeleton(nn).Dxyz(1,ff)+60; 
    end
    for nn = 1:length(compare_with_opt_skel) 
        not_optimized_skel(nn).Dxyz(2,ff)=not_optimized_skel(nn).Dxyz(2,ff)+15; 
        not_optimized_skel(nn).Dxyz(1,ff)=not_optimized_skel(nn).Dxyz(1,ff)+30;
    end
    for nn = 1:length(compare_with_opt_skel) 
        compare_with_not_opt_skel(nn).Dxyz(1,ff)=compare_with_not_opt_skel(nn).Dxyz(1,ff)+30; 
    end
    
%%%%%%%%%%
%RETARGET%
%%%%%%%%%%
    %%
    %left leg
    p1=[optimized_skel(11).Dxyz(3,ff);optimized_skel(11).Dxyz(1,ff);optimized_skel(11).Dxyz(2,ff)]; 
    p2=[optimized_skel(10).Dxyz(3,ff);optimized_skel(10).Dxyz(1,ff);optimized_skel(10).Dxyz(2,ff)]; 
    p3=[optimized_skel(9).Dxyz(3,ff);optimized_skel(9).Dxyz(1,ff);optimized_skel(9).Dxyz(2,ff)]; 
    all_angle=[0 0;];
    all_pt=[compare_with_opt_skel(11).Dxyz(3,ff) compare_with_opt_skel(11).Dxyz(1,ff) compare_with_opt_skel(11).Dxyz(2,ff);];
    displace_11_to_12=[optimized_skel(12).Dxyz(1,ff)-optimized_skel(11).Dxyz(1,ff) optimized_skel(12).Dxyz(2,ff)-optimized_skel(11).Dxyz(2,ff) optimized_skel(12).Dxyz(3,ff)-optimized_skel(11).Dxyz(3,ff)];
    displace_12_to_13=[optimized_skel(13).Dxyz(1,ff)-optimized_skel(12).Dxyz(1,ff) optimized_skel(13).Dxyz(2,ff)-optimized_skel(12).Dxyz(2,ff) optimized_skel(13).Dxyz(3,ff)-optimized_skel(12).Dxyz(3,ff)];
    %options = optimset('MaxFunEvals',1000000000000000000000000,'MaxIter',100000000000000000000);
    % A=[]; b=[];
    % xXx = fmincon(@angel_find_function_4,all_angle,A,b);
    xXx = fminsearch(@angel_find_function_4_3debug,all_angle);

    for i=1:1
        pt=[0;0;0];
        pt(1) = all_pt(i,1);
        pt(2) = all_pt(i,2);
        pt(3) = all_pt(i,3);
        if (((p1(1)-p2(1)==0)&&(p2(1)-p3(1)==0)&&(p1(1)-p3(1)==0)&&(pt(1)-p1(1)==0))||((p1(2)-p2(2)==0)&&(p2(2)-p3(2)==0)&&(p1(2)-p3(2)==0)&&(pt(2)-p1(2)==0))||((p1(3)-p2(3)==0)&&(p2(3)-p3(3)==0)&&(p1(3)-p3(3)==0)&&(pt(3)-p1(3)==0)))
            pts = [ p1' 1; p2' ,1; p3' 1]';
            ptspi =pts;
        else
            p3v = p3 + [0; 0; 5];
        % now p3v-p3-p2 form the original plane


        normal=cross(p3v-p3,p1-p3);  normal2=cross(p3v-p3,pt-p3);
        up=abs(normal(1).*normal2(1)+normal(2).*normal2(2)+normal(3).*normal2(3));
        down=sqrt( normal(1).*normal(1)+normal(2).*normal(2)+normal(3).*normal(3))*sqrt( normal2(1).*normal2(1)+normal2(2).*normal2(2)+normal2(3).*normal2(3) );
        % we want to rotate about Z-axis first
        radi= acos (up/down);
        axi=[p3(1)-p2(1);p3(2)-p2(2);p3(3)-p2(3)];
        ax2=axi/sqrt(sum(axi.^2));

        % rotate about vertical axis at pt3

        if(pt(2)<p2(2))
            M = AxelRot((radi) / pi * 180, [0; 0; 1], p3);
        elseif(pt(2)>p2(2))
            M = AxelRot((radi) / pi * 180, [0; 0; -1], p3);
        elseif(p3(2)<p2(2))
            M = AxelRot((radi) / pi * 180, [0; 0; -1], p3);
        else
            M = AxelRot((radi) / pi * 180, [0; 0; 1], p3);
        end

        Mx=sum(M([1 2 3],1)); 
        My=sum(M([1 2 3],2));
        Mz=sum(M([1 2 3],3));

        % construct points matrix
        pts = [ p1' 1; p2' ,1; p3' 1]';
        % transform points with M
        ptspi = M * pts;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        p1x=ptspi(1,1);p1y=ptspi(2,1); p1z=ptspi(3,1);p11=[p1x;p1y;p1z];
        p2x=ptspi(1,2);p2y=ptspi(2,2); p2z=ptspi(3,2);p22=[p2x;p2y;p2z];
        p3x=ptspi(1,3);p3y=ptspi(2,3); p3z=ptspi(3,3);p33=[p3x;p3y;p3z];
        normp=cross(p33-p22,p11-p22);%plane equation
        uninormp =normp/norm(normp);%unit vector

        ppp1=AxelRot((xXx(i,1))/ pi * 180, uninormp, p22);
        pt1pii = ppp1*[p11;1];

        % now rotate first joint (controlling p2-p1)
        xformPts = [ pt1pii(1:3) p22 p33];
        % now rotate wholr joint (controlling p3-p1)
        ppp1=AxelRot((xXx(i,2)) / pi * 180, uninormp, p33);
        ptspii = ppp1 * [xformPts [1; 1; 1]; 1 1 1 1 ];
        optimized_skel(11).Dxyz(3,ff)=ptspii(1,1);optimized_skel(11).Dxyz(1,ff)=ptspii(2,1);optimized_skel(11).Dxyz(2,ff)=ptspii(3,1);     
        optimized_skel(10).Dxyz(3,ff)=ptspii(1,2);optimized_skel(10).Dxyz(1,ff)=ptspii(2,2);optimized_skel(10).Dxyz(2,ff)=ptspii(3,2);
        optimized_skel(9).Dxyz(3,ff)=ptspii(1,3); optimized_skel(9).Dxyz(1,ff)=ptspii(2,3); optimized_skel(9).Dxyz(2,ff)=ptspii(3,3);
        optimized_skel(12).Dxyz(1,ff)=optimized_skel(11).Dxyz(1,ff)+displace_11_to_12(1);optimized_skel(12).Dxyz(2,ff)=optimized_skel(11).Dxyz(2,ff)+displace_11_to_12(2);optimized_skel(12).Dxyz(3,ff)=optimized_skel(11).Dxyz(3,ff)+displace_11_to_12(3);
        optimized_skel(13).Dxyz(1,ff)=optimized_skel(12).Dxyz(1,ff)+displace_12_to_13(1);optimized_skel(13).Dxyz(2,ff)=optimized_skel(12).Dxyz(2,ff)+displace_12_to_13(2);optimized_skel(13).Dxyz(3,ff)=optimized_skel(12).Dxyz(3,ff)+displace_12_to_13(3);
    end

    %%
    %right leg
    
    p1=[optimized_skel(5).Dxyz(3,ff);optimized_skel(5).Dxyz(1,ff);optimized_skel(5).Dxyz(2,ff)]; 
    p2=[optimized_skel(4).Dxyz(3,ff);optimized_skel(4).Dxyz(1,ff);optimized_skel(4).Dxyz(2,ff)]; 
    p3=[optimized_skel(3).Dxyz(3,ff);optimized_skel(3).Dxyz(1,ff);optimized_skel(3).Dxyz(2,ff)]; 
    displace_5_to_6=[optimized_skel(6).Dxyz(1,ff)-optimized_skel(5).Dxyz(1,ff) optimized_skel(6).Dxyz(2,ff)-optimized_skel(5).Dxyz(2,ff) optimized_skel(6).Dxyz(3,ff)-optimized_skel(5).Dxyz(3,ff)];
    displace_6_to_7=[optimized_skel(7).Dxyz(1,ff)-optimized_skel(6).Dxyz(1,ff) optimized_skel(7).Dxyz(2,ff)-optimized_skel(6).Dxyz(2,ff) optimized_skel(7).Dxyz(3,ff)-optimized_skel(6).Dxyz(3,ff)];
    all_pt=[compare_with_opt_skel(5).Dxyz(3,ff) compare_with_opt_skel(5).Dxyz(1,ff) compare_with_opt_skel(5).Dxyz(2,ff);];
    xXx = fminsearch(@angel_find_function_4_3debug,all_angle);
    % A=[]; b=[];
    % xXx = fmincon(@angel_find_function_4,all_angle,A,b);
    
    for i=1:1
    pt=[0;0;0];
    pt(1) = all_pt(i,1);
    pt(2) = all_pt(i,2);
    pt(3) = all_pt(i,3);

    if (((p1(1)-p2(1)==0)&&(p2(1)-p3(1)==0)&&(p1(1)-p3(1)==0)&&(pt(1)-p1(1)==0))||((p1(2)-p2(2)==0)&&(p2(2)-p3(2)==0)&&(p1(2)-p3(2)==0)&&(pt(2)-p1(2)==0))||((p1(3)-p2(3)==0)&&(p2(3)-p3(3)==0)&&(p1(3)-p3(3)==0)&&(pt(3)-p1(3)==0)))
        pts = [ p1' 1; p2' ,1; p3' 1]';
        ptspi =pts;
    else
        p3v = p3 + [0; 0; 5];
    % now p3v-p3-p2 form the original plane

    normal=cross(p3v-p3,p1-p3);  normal2=cross(p3v-p3,pt-p3);
    up=abs(normal(1).*normal2(1)+normal(2).*normal2(2)+normal(3).*normal2(3)  );
    down=sqrt(normal(1).*normal(1)+normal(2).*normal(2)+normal(3).*normal(3))*sqrt(normal2(1).*normal2(1)+normal2(2).*normal2(2)+normal2(3).*normal2(3));
    % we want to rotate about Z-axis first
    radi= acos (up/down);
    %axi=[normal(1);normal(2);normal(3)];
    axi=[p3(1)-p2(1);p3(2)-p2(2);p3(3)-p2(3)];
    ax2=axi/sqrt(sum(axi.^2));
    % rotate about vertical axis at pt3

    if(pt(2)<p2(2))
        M = AxelRot((radi) / pi * 180, [0; 0; 1], p3);
    elseif(pt(2)>p2(2))
        M = AxelRot((radi) / pi * 180, [0; 0; -1], p3);
    elseif(p3(2)<p2(2))
        M = AxelRot((radi) / pi * 180, [0; 0; -1], p3);
    else
        M = AxelRot((radi) / pi * 180, [0; 0; 1], p3);
    end

    Mx=sum(M([1 2 3],1)); 
    My=sum(M([1 2 3],2));
    Mz=sum(M([1 2 3],3));

    % construct points matrix
    pts = [ p1' 1; p2' ,1; p3' 1]';
    ptspi = M * pts;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p1x=ptspi(1,1);p1y=ptspi(2,1); p1z=ptspi(3,1);p11=[p1x;p1y;p1z];
    p2x=ptspi(1,2);p2y=ptspi(2,2); p2z=ptspi(3,2);p22=[p2x;p2y;p2z];
    p3x=ptspi(1,3);p3y=ptspi(2,3); p3z=ptspi(3,3);p33=[p3x;p3y;p3z];
    normp=cross(p33-p22,p11-p22);
    uninormp =normp/norm(normp);

    ppp1=AxelRot((xXx(i,1))/ pi * 180, uninormp, p22);
    pt1pii = ppp1*[p11;1];

    % now rotate first joint (controlling p2-p1)
    xformPts = [ pt1pii(1:3) p22 p33];
    % now rotate wholr joint (controlling p3-p1)
    ppp1=AxelRot((xXx(i,2)) / pi * 180, uninormp, p33);
    ptspii = ppp1 * [xformPts [1; 1; 1]; 1 1 1 1 ];
    optimized_skel(5).Dxyz(3,ff)=ptspii(1,1);optimized_skel(5).Dxyz(1,ff)=ptspii(2,1);optimized_skel(5).Dxyz(2,ff)=ptspii(3,1);     
    optimized_skel(4).Dxyz(3,ff)=ptspii(1,2);optimized_skel(4).Dxyz(1,ff)=ptspii(2,2);optimized_skel(4).Dxyz(2,ff)=ptspii(3,2);
    optimized_skel(3).Dxyz(3,ff)=ptspii(1,3);optimized_skel(3).Dxyz(1,ff)=ptspii(2,3);optimized_skel(3).Dxyz(2,ff)=ptspii(3,3);
    optimized_skel(6).Dxyz(1,ff)=optimized_skel(5).Dxyz(1,ff)+displace_5_to_6(1);optimized_skel(6).Dxyz(2,ff)=optimized_skel(5).Dxyz(2,ff)+displace_5_to_6(2);optimized_skel(6).Dxyz(3,ff)=optimized_skel(5).Dxyz(3,ff)+displace_5_to_6(3);
    optimized_skel(7).Dxyz(1,ff)=optimized_skel(6).Dxyz(1,ff)+displace_6_to_7(1);optimized_skel(7).Dxyz(2,ff)=optimized_skel(6).Dxyz(2,ff)+displace_6_to_7(2);optimized_skel(7).Dxyz(3,ff)=optimized_skel(6).Dxyz(3,ff)+displace_6_to_7(3);
    end

    %%
    for x=-40:40:40 
        for y=-40:40:40
            plot3([x x],[y y] ,[-40 40],'r-','LineWidth',0.1);
        end 
    end
    for y=-40:40:40 
        for z=-40:40:40
            plot3([-40 40],[y y] ,[z z],'r-','LineWidth',0.1);
        end
    end
    for z=-40:40:40 
        for x=-40:40:40
            plot3([x x],[-40 40] ,[z z],'r-','LineWidth',0.1);
        end
    end

%%
    for nn = 1:length(compare_with_opt_skel) 
        left_colour_change_begin=5;
        right_colour_change_begin=11;
        hold on
        
        compare_opt_skel_parent = compare_with_opt_skel(nn).compare_opt_skel_parent;
        opt_skel_parent = optimized_skel(nn).compare_opt_skel_parent;
        original_skel_parent = original_skeleton(nn).compare_opt_skel_parent;
        not_opt_skel_parent = not_optimized_skel(nn).compare_opt_skel_parent;
        compare_not_opt_skel_skeleton = compare_with_not_opt_skel(nn).compare_opt_skel_parent;
        
        if((nn==11)||(nn==12)||(nn==13)||(nn==5)||(nn==6)||(nn==7))
            plot3(compare_with_opt_skel(nn).Dxyz(1,ff),compare_with_opt_skel(nn).Dxyz(3,ff),compare_with_opt_skel(nn).Dxyz(2,ff),'g.','markersize',10)
            plot3(optimized_skel(nn).Dxyz(1,ff),optimized_skel(nn).Dxyz(3,ff),optimized_skel(nn).Dxyz(2,ff),'g.','markersize',10)
            plot3(original_skeleton(nn).Dxyz(1,ff),original_skeleton(nn).Dxyz(3,ff),original_skeleton(nn).Dxyz(2,ff),'g.','markersize',10)
            plot3(not_optimized_skel(nn).Dxyz(1,ff),not_optimized_skel(nn).Dxyz(3,ff),not_optimized_skel(nn).Dxyz(2,ff),'g.','markersize',10)
            plot3(compare_with_not_opt_skel(nn).Dxyz(1,ff),compare_with_not_opt_skel(nn).Dxyz(3,ff),compare_with_not_opt_skel(nn).Dxyz(2,ff),'g.','markersize',10)
        else
            plot3(compare_with_opt_skel(nn).Dxyz(1,ff),compare_with_opt_skel(nn).Dxyz(3,ff),compare_with_opt_skel(nn).Dxyz(2,ff),'y.','markersize',20)
            plot3(optimized_skel(nn).Dxyz(1,ff),optimized_skel(nn).Dxyz(3,ff),optimized_skel(nn).Dxyz(2,ff),'y.','markersize',20)
            plot3(original_skeleton(nn).Dxyz(1,ff),original_skeleton(nn).Dxyz(3,ff),original_skeleton(nn).Dxyz(2,ff),'y.','markersize',20)
            plot3(not_optimized_skel(nn).Dxyz(1,ff),not_optimized_skel(nn).Dxyz(3,ff),not_optimized_skel(nn).Dxyz(2,ff),'y.','markersize',20)
            plot3(compare_with_not_opt_skel(nn).Dxyz(1,ff),compare_with_not_opt_skel(nn).Dxyz(3,ff),compare_with_not_opt_skel(nn).Dxyz(2,ff),'y.','markersize',20)
        end
        
        grid minor

        %lines, dxyz is point position
        if compare_opt_skel_parent > 0
            plot3([compare_with_opt_skel(compare_opt_skel_parent).Dxyz(1,ff) compare_with_opt_skel(nn).Dxyz(1,ff)],...
            [compare_with_opt_skel(compare_opt_skel_parent).Dxyz(3,ff) compare_with_opt_skel(nn).Dxyz(3,ff)],...
            [compare_with_opt_skel(compare_opt_skel_parent).Dxyz(2,ff) compare_with_opt_skel(nn).Dxyz(2,ff)],'k-','LineWidth',3) 
        end

        plot3([compare_with_opt_skel(right_colour_change_begin).Dxyz(1,ff) compare_with_opt_skel(13).Dxyz(1,ff)],...
        [compare_with_opt_skel(right_colour_change_begin).Dxyz(3,ff) compare_with_opt_skel(13).Dxyz(3,ff)],...
        [compare_with_opt_skel(right_colour_change_begin).Dxyz(2,ff) compare_with_opt_skel(13).Dxyz(2,ff)],'r-','LineWidth',3)
        plot3([compare_with_opt_skel(left_colour_change_begin).Dxyz(1,ff) compare_with_opt_skel(7).Dxyz(1,ff)],...
        [compare_with_opt_skel(left_colour_change_begin).Dxyz(3,ff) compare_with_opt_skel(7).Dxyz(3,ff)],...
        [compare_with_opt_skel(left_colour_change_begin).Dxyz(2,ff) compare_with_opt_skel(7).Dxyz(2,ff)],'r-','LineWidth',3)

        if opt_skel_parent > 0
            plot3([optimized_skel(opt_skel_parent).Dxyz(1,ff) optimized_skel(nn).Dxyz(1,ff)],...
            [optimized_skel(opt_skel_parent).Dxyz(3,ff) optimized_skel(nn).Dxyz(3,ff)],...
            [optimized_skel(opt_skel_parent).Dxyz(2,ff) optimized_skel(nn).Dxyz(2,ff)],'g-','LineWidth',3)
            grid minor
        end
        if original_skel_parent > 0
            plot3([original_skeleton(original_skel_parent).Dxyz(1,ff) original_skeleton(nn).Dxyz(1,ff)],...
            [original_skeleton(original_skel_parent).Dxyz(3,ff) original_skeleton(nn).Dxyz(3,ff)],...
            [original_skeleton(original_skel_parent).Dxyz(2,ff) original_skeleton(nn).Dxyz(2,ff)],'k-','LineWidth',3)
            grid minor
        end
        plot3([original_skeleton(right_colour_change_begin).Dxyz(1,ff) original_skeleton(13).Dxyz(1,ff)],...
        [original_skeleton(right_colour_change_begin).Dxyz(3,ff) original_skeleton(13).Dxyz(3,ff)],...
        [original_skeleton(right_colour_change_begin).Dxyz(2,ff) original_skeleton(13).Dxyz(2,ff)],'r-','LineWidth',3)
        plot3([original_skeleton(left_colour_change_begin).Dxyz(1,ff) original_skeleton(7).Dxyz(1,ff)],...
        [original_skeleton(left_colour_change_begin).Dxyz(3,ff) original_skeleton(7).Dxyz(3,ff)],...
        [original_skeleton(left_colour_change_begin).Dxyz(2,ff) original_skeleton(7).Dxyz(2,ff)],'r-','LineWidth',3)

        if not_opt_skel_parent > 0
            plot3([not_optimized_skel(not_opt_skel_parent).Dxyz(1,ff) not_optimized_skel(nn).Dxyz(1,ff)],...
            [not_optimized_skel(not_opt_skel_parent).Dxyz(3,ff) not_optimized_skel(nn).Dxyz(3,ff)],...
            [not_optimized_skel(not_opt_skel_parent).Dxyz(2,ff) not_optimized_skel(nn).Dxyz(2,ff)],'g-','LineWidth',3)
            grid minor
        end

        if compare_not_opt_skel_skeleton > 0
            plot3([compare_with_not_opt_skel(compare_not_opt_skel_skeleton).Dxyz(1,ff) compare_with_not_opt_skel(nn).Dxyz(1,ff)],...
            [compare_with_not_opt_skel(compare_not_opt_skel_skeleton).Dxyz(3,ff) compare_with_not_opt_skel(nn).Dxyz(3,ff)],...
            [compare_with_not_opt_skel(compare_not_opt_skel_skeleton).Dxyz(2,ff) compare_with_not_opt_skel(nn).Dxyz(2,ff)],'k-','LineWidth',3)
            grid minor
        end

        plot3([compare_with_not_opt_skel(right_colour_change_begin).Dxyz(1,ff) compare_with_not_opt_skel(13).Dxyz(1,ff)],...
        [compare_with_not_opt_skel(right_colour_change_begin).Dxyz(3,ff) compare_with_not_opt_skel(13).Dxyz(3,ff)],...
        [compare_with_not_opt_skel(right_colour_change_begin).Dxyz(2,ff) compare_with_not_opt_skel(13).Dxyz(2,ff)],'r-','LineWidth',3)
        plot3([compare_with_not_opt_skel(left_colour_change_begin).Dxyz(1,ff) compare_with_not_opt_skel(7).Dxyz(1,ff)],...
        [compare_with_not_opt_skel(left_colour_change_begin).Dxyz(3,ff) compare_with_not_opt_skel(7).Dxyz(3,ff)],...
        [compare_with_not_opt_skel(left_colour_change_begin).Dxyz(2,ff) compare_with_not_opt_skel(7).Dxyz(2,ff)],'r-','LineWidth',3)
        
    end
    
    view(-135,30)
    axis equal on 

    drawnow
    if write_video, writeVideo(vidObj,getframe); end

end

if write_video, close(vidObj); end