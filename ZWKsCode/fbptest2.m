% [fn,pn,fi]=uigetfile('*.bmp;*.jpg','选择图片');
% I=imread([pn fn]);
% P1=I(:,:,1);
 P1=phantom( 128 );
sp=[200,200];%first position of projection source
Sr0=sp(1);
Sc0=sp(2);
[r,c]=size(P1);
P1=double(P1);
newP1=zeros(r,c);
L=ceil(sqrt(r^2+c^2));%long axis of imaging
%y=(r/c)*(x-r/2)+c/2; %中垂线方程
%x=solve('x^2+(64+(x-64))^2=200^2*2','x');
r_mov=(r-1)/2;
c_mov=(c-1)/2;%将图像中心移到坐标原点
d_mov=(r+1)/2;%移动的距离
theta=acosd(dot([r_mov-Sr0,-c_mov-Sc0],[-r_mov-Sr0,c_mov-Sc0])/(norm([r_mov-Sr0,-c_mov-Sc0])*norm([-r_mov-Sr0,c_mov-Sc0])));%minimum angle of fan beam
theta=ceil(theta*1.2);%此时theta为51
proj_num=220;%每个扇束分128个角度,此时射线有129条,包括中心线和-64,64的两侧线
p_theta=theta/proj_num;
filterR_L=zeros(proj_num+1,1);
rot_num=359;%转动次数127,再加原位置1次共128次
R=zeros(proj_num+1,rot_num+1);
R_VIEW=R;
rot_angle=1;%每次转动射线源角度
%%different position of projection source
for t=0:rot_num %投影源转动rot_num次，每次rot_angle度，共180+theta度
    Sr=cos(rot_angle*t*pi/180)*Sr0-sin(rot_angle*t*pi/180)*Sc0;
    Sc=sin(rot_angle*t*pi/180)*Sr0+cos(rot_angle*t*pi/180)*Sc0;
    %     Sr=cos(180*pi/180)*Sr0-sin(180*pi/180)*Sc0;
    %     Sc=sin(180*pi/180)*Sr0+cos(180*pi/180)*Sc0;
    % Sr=-4;
    %Sc=-4;
    s=0;
%     if t==64
%         pause1=1;
%     end
    for n=-proj_num/2:proj_num/2
        %n=-35;
        %gamma=45+180+theta/2-(n+75)*p_theta;
        gamma=45+rot_angle*t+theta/2-(n+proj_num/2+1)*p_theta;%其中，中心线是n=0的投影线。顺时针转动射线，逆时针转动射线源
        gamma=gamma*pi/180;%转换成弧度
        s=s+1;%s=1是第一根射线
        k=0;
%         if s==54
%             pause2=1;
%         end
        %直线tan(gamma)*(x-Sr)+Sc-y=0
        cross1=roundn((Sc+tan(gamma)*(-(r-1)/2-Sr)+(c-1)/2)*(Sc+tan(gamma)*((r-1)/2-Sr)+(c-1)/2),-4);%线段1端点(-(r-1)/2,-(c-1)/2)和((r-1)/2,-(c-1)/2)
        cross2=roundn((Sc+tan(gamma)*((r-1)/2-Sr)+(c-1)/2)*(Sc+tan(gamma)*((r-1)/2-Sr)-(c-1)/2),-4);%线段2端点((r-1)/2,-(c-1)/2)和((r-1)/2,(c-1)/2)
        cross3=roundn((Sc+tan(gamma)*((r-1)/2-Sr)-(c-1)/2)*(Sc+tan(gamma)*(-(r-1)/2-Sr)-(c-1)/2),-4);%线段3端点((r-1)/2,(c-1)/2)和(-(r-1)/2,(c-1)/2)
        cross4=roundn((Sc+tan(gamma)*(-(r-1)/2-Sr)-(c-1)/2)*(Sc+tan(gamma)*(-(r-1)/2-Sr)+(c-1)/2),-4);%线段4端点(-(r-1)/2,(c-1)/2)和(-(r-1)/2,-(c-1)/2)
        if  tan(gamma)>10^8||tan(gamma)<-10^8
            if Sr>-r_mov&&Sr<r_mov
                Sr_fore=floor(Sr+d_mov);
                for y=1:c-1
                    k=k+P1(Sr_fore,y);
                end
            end
        elseif -10^-5<tan(gamma)&&tan(gamma)<10^-5
            if Sc>-r_mov&&Sc<r_mov
                Sc_fore=floor(roundn(Sc,-4)+d_mov);
                for x=1:r-1
                    k=k+P1(x,Sc_fore);
                end
            end
        elseif tan(gamma)>0
            
            %%第一种情况，射线过下边
            
            if cross1<=0
                xp1=(-(c-1)/2-Sc)/tan(gamma)+Sr;%与下边的交点坐标(xp1,-(c-1)/2)
                xp1=roundn(xp1,-4);
                xp1_movback=xp1+d_mov;%映射回原图像
                xp1_fore=floor(xp1_movback);%点（xp1,-(c-1)/2）所对应的像素在原图中为（xp1_fore，1），作为起始像素点
                %射线过下边的同时过右边，并且不能同时过右下角同一个点((r-1)/2,-(c-1)/2)，即(r_mov,-c_mov)
                if cross2<=0&&(-c_mov~=tan(gamma)*(r_mov-Sc)+Sr)
                    xp=zeros(r-xp1_fore+1,1);%射线与像素网格相交的点数
                    yp=zeros(r-xp1_fore+1,1);
                    for i=2:r-xp1_fore+1
                        xp(i)=xp1_fore+i-1-d_mov;%第i个点的横坐标
                        yp(i)=Sc+tan(gamma)*(xp(i)-Sr);%第i个点的纵坐标
                    end
                    xp(1)=xp1;
                    yp(1)=-c_mov;
                    for i=1:r-xp1_fore
                        xp_fore=floor(roundn(xp(i),-4)+d_mov)-d_mov;
                        yp_fore=floor(roundn(yp(i),-4)+d_mov)-d_mov;
                        yp_ceil=ceil(roundn(yp(i+1),-4)+d_mov)-d_mov;
                        dot_num=yp_ceil-yp_fore;%相邻两点之间穿过的像素个数，也等于穿过的射线分段数
                        if dot_num>1
                            x_dot=zeros(dot_num-1,1);
                            y_dot=zeros(dot_num-1,1);
                            for j=1:dot_num-1
                                y_dot(j)=yp_fore+j;%xp(i)与xp(i+1)之间的穿过的像素的点的纵坐标
                                x_dot(j)=(y_dot(j)-Sc)/tan(gamma)+Sr;%xp(i)与xp(i+1)之间的穿过的像素的点的横坐标
                            end
                            d=zeros(dot_num,1);
                            if dot_num>2
                                for d_count=2:dot_num-1
                                    d(1)=sqrt((xp(i)-x_dot(1))^2+(yp(i)-y_dot(1))^2);%第i个点与第i+1个点之间穿过像素第一段的距离
                                    d(d_count)=sqrt((x_dot(d_count-1)-x_dot(d_count))^2+(y_dot(d_count-1)-y_dot(d_count))^2);
                                    d(dot_num)=sqrt((x_dot(dot_num-1)-xp(i+1))^2+(y_dot(dot_num-1)-yp(i+1))^2);
                                end
                            else
                                d(1)=sqrt((xp(i)-x_dot(1))^2+(yp(i)-y_dot(1))^2);
                                d(2)=sqrt((x_dot(dot_num-1)-xp(i+1))^2+(y_dot(dot_num-1)-yp(i+1))^2);
                            end
                            for k_count=1:dot_num
                                k=k+P1(xp_fore+d_mov,yp_fore+d_mov+k_count-1)*d(k_count);
                            end
                        else
                            d0=sqrt((xp(i)-xp(i+1))^2+(yp(i)-yp(i+1))^2);
                            k=k+P1(xp_fore+d_mov,yp_fore+d_mov)*d0;
                        end
                    end
                end
                
                %射线过下边的同时过上边
                if cross3<0
                    x_end=(c_mov-Sc)/tan(gamma)+Sr;%与上边界交点(x_end,c_mov)
                    x_endmovback=x_end+d_mov;
                    x_endfloor=floor(x_endmovback);
                    if x_endfloor==x_endmovback %最后一个点横坐标是整数
                        num=x_endfloor-xp1_fore+1;
                        for i=2:num
                            xp(i)=xp1_fore+i-1-d_mov;%第i个点的横坐标
                            yp(i)=Sc+tan(gamma)*(xp(i)-Sr);%第i个点的纵坐标
                        end
                        xp(1)=xp1;
                        yp(1)=-c_mov;
                    else %最后一个点横坐标不是整数
                        num=x_endfloor-xp1_fore+2;
                        if num<3
                            xp(1)=xp1;
                            xp(2)=x_end;
                            yp(1)=-c_mov;
                            yp(2)=c_mov;
                        else
                            for i=2:num-1
                                xp(i)=xp1_fore+i-1-d_mov;%第i个点的横坐标
                                yp(i)=Sc+tan(gamma)*(xp(i)-Sr);%第i个点的纵坐标
                            end
                            xp(1)=xp1;
                            yp(1)=-c_mov;
                            xp(num)=x_end;
                            yp(num)=c_mov;
                        end
                    end
                    for i=1:num-1
                        xp_fore=floor(roundn(xp(i),-4)+d_mov)-d_mov;
                        yp_fore=floor(roundn(yp(i),-4)+d_mov)-d_mov;
                        yp_ceil=ceil(roundn(yp(i+1),-4)+d_mov)-d_mov;
                        dot_num=yp_ceil-yp_fore;%相邻两点之间穿过的像素个数，也等于穿过的射线分段数
                        if dot_num>1
                            x_dot=zeros(dot_num-1,1);
                            y_dot=zeros(dot_num-1,1);
                            for j=1:dot_num-1
                                y_dot(j)=yp_fore+j;%xp(i)与xp(i+1)之间的穿过的像素的点的纵坐标
                                x_dot(j)=(y_dot(j)-Sc)/tan(gamma)+Sr;%xp(i)与xp(i+1)之间的穿过的像素的点的横坐标
                            end
                            d=zeros(dot_num,1);
                            if dot_num>2
                                for d_count=2:dot_num-1
                                    d(1)=sqrt((xp(i)-x_dot(1))^2+(yp(i)-y_dot(1))^2);%第i个点与第i+1个点之间穿过像素第一段的距离
                                    d(d_count)=sqrt((x_dot(d_count-1)-x_dot(d_count))^2+(y_dot(d_count-1)-y_dot(d_count))^2);
                                    d(dot_num)=sqrt((x_dot(dot_num-1)-xp(i+1))^2+(y_dot(dot_num-1)-yp(i+1))^2);
                                end
                            else
                                d(1)=sqrt((xp(i)-x_dot(1))^2+(yp(i)-y_dot(1))^2);
                                d(dot_num)=sqrt((x_dot(dot_num-1)-xp(i+1))^2+(y_dot(dot_num-1)-yp(i+1))^2);
                            end
                            for k_count=1:dot_num
                                if xp_fore==xp;
                                    xp_fore=xp_fore-1;
                                end
                                k=k+P1(xp_fore+d_mov,yp_fore+k_count-1+d_mov)*d(k_count);
                            end
                        else
                            if xp_fore==xp;
                                xp_fore=xp_fore-1;
                            end
                            d0=sqrt((xp(i)-xp(i+1))^2+(yp(i)-yp(i+1))^2);
                            k=k+P1(xp_fore+d_mov,yp_fore+d_mov)*d0;
                        end
                    end
                    
                end
            
            
            %%第二种情况，射线过左边
            
        elseif cross4<0
                if cross2<=0 %与左边相交的情况下与右边相交
                    xp4=zeros(r,1);
                    yp4=zeros(r,1);
                    for i=1:r
                        xp4(i)=i-d_mov;
                        yp4(i)=tan(gamma)*(xp4(i)-Sr)+Sc;%射线与x=i相交的交点(xp4,yp4)
                    end
                    for i=1:r-1
                        xp4_fore=floor(roundn(xp4(i),-4)+d_mov)-d_mov;
                        yp4_fore=floor(roundn(yp4(i),-4)+d_mov)-d_mov;
                        yp4_back=floor(roundn(yp4(i+1),-4)+d_mov)-d_mov;
                        if yp4_fore==yp4_back %说明相邻两个点在同一个像素区域内
                            d4=sqrt(1+(yp4(i)-yp4(i+1))^2);
                            k=k+P1(xp4_fore+d_mov,yp4_fore+d_mov)*d4;
                        else %说明相邻两个点在相邻两个像素区域内
                            yp4_mid=floor(roundn(yp4(i+1),-4)+d_mov)-d_mov;
                            xp4_mid=(yp4_mid-Sc)/tan(gamma)+Sr;
                            d41=sqrt((xp4(i)-xp4_mid)^2+(yp4(i)- yp4_mid)^2);
                            d42=sqrt((xp4(i+1)-xp4_mid)^2+(yp4(i+1)- yp4_mid)^2);
                            k=k+P1(xp4_fore+d_mov,yp4_fore+d_mov)*d41+P1(xp4_fore+d_mov,yp4_mid+d_mov)*d42;
                        end
                    end
                end
                if cross3<0 %与左边相交的情况下与上边相交
                    xp4_end=(c_mov-Sc)/tan(gamma)+Sr; %射线与上边交点(xp4,c_mov)
                    num=ceil(roundn(xp4_end,-4)+d_mov); %射线与图像纵轴线相交的点数
                    xp4=zeros(num,1);
                    yp4=zeros(num,1);
                    %获取与图像纵轴线相交的点所有坐标
                    for i=1:num-1
                        xp4(i)=i-d_mov;
                        yp4(i)=tan(gamma)*(xp4(i)-Sr)+Sc;
                    end
                    xp4(num)=xp4_end;
                    yp4(num)=c_mov;
                    for i=1:num-1
                        xp4_fore=floor(roundn(xp4(i),-4)+d_mov)-d_mov;
                        yp4_fore=floor(roundn(yp4(i),-4)+d_mov)-d_mov;
                        yp4_ceil=ceil(roundn(yp4(i+1),-4)+d_mov)-d_mov;
                        dot_num=yp4_ceil-yp4_fore; %与图像纵轴线相交的相邻两个点之间穿过的像素个数，也等于穿过的射线分段数
                        x_dot=zeros(dot_num-1,1);
                        y_dot=zeros(dot_num-1,1);
                        d=zeros(dot_num,1);
                        if dot_num>1
                            for j=1:dot_num-1
                                y_dot(j)=yp4_fore+j;%xp(i)与xp(i+1)之间的穿过的像素的点的纵坐标
                                x_dot(j)=(y_dot(j)-Sc)/tan(gamma)+Sr;%xp(i)与xp(i+1)之间的穿过的像素的点的横坐标
                            end
                            d=zeros(dot_num,1);
                            if dot_num>2
                                for d_count=2:dot_num-1
                                    d(1)=sqrt((xp4(i)-x_dot(1))^2+(yp4(i)-y_dot(1))^2);%第i个点与第i+1个点之间穿过像素第一段的距离
                                    d(d_count)=sqrt((x_dot(d_count-1)-x_dot(d_count))^2+(y_dot(d_count-1)-y_dot(d_count))^2);
                                    d(dot_num)=sqrt((x_dot(dot_num-1)-xp4(i+1))^2+(y_dot(dot_num-1)-yp4(i+1))^2);
                                end
                            else
                                d(1)=sqrt((xp4(i)-x_dot(1))^2+(yp4(i)-y_dot(1))^2);
                                d(dot_num)=sqrt((x_dot(dot_num-1)-xp4(i+1))^2+(y_dot(dot_num-1)-yp4(i+1))^2);
                            end
                        else
                            d(1)=sqrt((xp4(i)-xp4(i+1))^2+(yp4(i)-yp4(i+1))^2);
                        end
                        for k_count=1:dot_num
                            k=k+P1(xp4_fore+d_mov,yp4_fore+k_count-1+d_mov)*d(k_count);
                        end
                    end
                end
            end
        
            
            
            
            %tan(gamma)<0情况
        elseif tan(gamma)<0 %%第三种情况，射线过下边
            if cross1<=0
                xp1=(-(c-1)/2-Sc)/tan(gamma)+Sr;;%与下边的交点坐标(xp1,-c_mov)
                xp1=roundn(xp1,-4);
                xp1_ceil=ceil(roundn(xp1,-4)+d_mov);%点（xp1,1）所在的像素为（xp1_fore，1），作为起始像素点
                xp=zeros(xp1_ceil,1);%射线与像素网格相交的点数是xp1_fore
                yp=zeros(xp1_ceil,1);
                %射线过下边的同时过左边，并且不能同时过左下角同一个点(-r_mov,-c_mov)
                if cross4<=0&&(-c_mov~=tan(gamma)*(-r_mov-Sc)+Sr)
                    for i=2:xp1_ceil
                        xp(i)=xp1_ceil-i+1-d_mov;%第i个点的横坐标
                        yp(i)=Sc+tan(gamma)*(xp(i)-Sr);%第i个点的纵坐标
                    end
                    xp(1)=xp1;
                    yp(1)=-c_mov;
                    for i=1:xp1_ceil-1
                        xp_fore=floor(roundn(xp(i),-4)+d_mov)-d_mov;
                        yp_fore=floor(roundn(yp(i),-4)+d_mov)-d_mov;
                        yp_ceil=ceil(roundn(yp(i+1),-4)+d_mov)-d_mov;
                        dot_num=yp_ceil-yp_fore;%相邻两点之间穿过的像素个数，也等于穿过的射线分段数
                        if dot_num>1
                            x_dot=zeros(dot_num-1,1);
                            y_dot=zeros(dot_num-1,1);
                            for j=1:dot_num-1
                                y_dot(j)=yp_fore+j;%xp(i)与xp(i+1)之间的穿过的像素的点的纵坐标
                                x_dot(j)=(y_dot(j)-Sc)/tan(gamma)+Sr;%xp(i)与xp(i+1)之间的穿过的像素的点的横坐标
                            end
                            d=zeros(dot_num,1);
                            if dot_num>2
                                for d_count=2:dot_num-1
                                    d(1)=sqrt((xp(i)-x_dot(1))^2+(yp(i)-y_dot(1))^2);%第i个点与第i+1个点之间穿过像素第一段的距离
                                    d(d_count)=sqrt((x_dot(d_count-1)-x_dot(d_count))^2+(y_dot(d_count-1)-y_dot(d_count))^2);
                                    d(dot_num)=sqrt((x_dot(dot_num-1)-xp(i+1))^2+(y_dot(dot_num-1)-yp(i+1))^2);
                                end
                            else
                                d(1)=sqrt((xp(i)-x_dot(1))^2+(yp(i)-y_dot(1))^2);
                                d(2)=sqrt((x_dot(dot_num-1)-xp(i+1))^2+(y_dot(dot_num-1)-yp(i+1))^2);
                            end
                            for k_count=1:dot_num
                                if xp_fore==xp;
                                    xp_fore=xp_fore-1;
                                end
                                k=k+P1(xp_fore+d_mov,yp_fore+k_count-1+d_mov)*d(k_count);
                            end
                        else
                            if xp_fore==xp;
                                xp_fore=xp_fore-1;
                            end
                            d0=sqrt((xp(i)-xp(i+1))^2+(yp(i)-yp(i+1))^2);
                            k=k+P1(xp_fore+d_mov,yp_fore+d_mov)*d0;
                        end
                    end
                end
                
                %射线过下边的同时过上边
                if cross3<0
                    x_end=(c_mov-Sc)/tan(gamma)+Sr;%与上边界交点(x_end,c_mov)
                    x_endmovback=x_end+d_mov;
                    x_endceil=ceil(x_endmovback);
                    if x_endceil==x_endmovback %最后一个点横坐标是整数
                        num=xp1_ceil-x_endceil+1;
                        for i=2:num
                            xp(i)=xp1_ceil-i+1-d_mov;%第i个点的横坐标
                            yp(i)=Sc+tan(gamma)*(xp(i)-Sr);%第i个点的纵坐标
                        end
                        xp(1)=xp1;
                        yp(1)=-c_mov;
                    else %最后一个点横坐标不是整数
                        num=xp1_ceil-x_endceil+2;
                        if num<3
                            xp(1)=xp1;
                            xp(2)=x_end;
                            yp(1)=-c_mov;
                            yp(2)=c_mov;
                        else
                            for i=2:num-1
                                xp(i)=xp1_ceil-i+1-d_mov;%第i个点的横坐标
                                yp(i)=Sc+tan(gamma)*(xp(i)-Sr);%第i个点的纵坐标
                            end
                            xp(1)=xp1;
                            yp(1)=-c_mov;
                            xp(num)=x_end;
                            yp(num)=c_mov;
                        end
                    end
                    for i=1:num-1
                        xp_fore=floor(roundn(xp(i),-4)+d_mov)-d_mov;
                        yp_fore=floor(roundn(yp(i),-4)+d_mov)-d_mov;
                        yp_ceil=ceil(roundn(yp(i+1),-4)+d_mov)-d_mov;
                        dot_num=yp_ceil-yp_fore;%相邻两点之间穿过的像素个数，也等于穿过的射线分段数
                        if dot_num>1
                            x_dot=zeros(dot_num-1,1);
                            y_dot=zeros(dot_num-1,1);
                            for j=1:dot_num-1
                                y_dot(j)=yp_fore+j;%xp(i)与xp(i+1)之间的穿过的像素的点的纵坐标
                                x_dot(j)=(y_dot(j)-Sc)/tan(gamma)+Sr;%xp(i)与xp(i+1)之间的穿过的像素的点的横坐标
                            end
                            d=zeros(dot_num,1);
                            if dot_num>2
                                for d_count=2:dot_num-1
                                    d(1)=sqrt((xp(i)-x_dot(1))^2+(yp(i)-y_dot(1))^2);%第i个点与第i+1个点之间穿过像素第一段的距离
                                    d(d_count)=sqrt((x_dot(d_count-1)-x_dot(d_count))^2+(y_dot(d_count-1)-y_dot(d_count))^2);
                                    d(dot_num)=sqrt((x_dot(dot_num-1)-xp(i+1))^2+(y_dot(dot_num-1)-yp(i+1))^2);
                                end
                            else
                                d(1)=sqrt((xp(i)-x_dot(1))^2+(yp(i)-y_dot(1))^2);
                                d(dot_num)=sqrt((x_dot(dot_num-1)-xp(i+1))^2+(y_dot(dot_num-1)-yp(i+1))^2);
                            end
                            for k_count=1:dot_num
                                if xp_fore==xp;
                                    xp_fore=xp_fore-1;
                                end
                                k=k+P1(xp_fore+d_mov,yp_fore+k_count-1+d_mov)*d(k_count);
                            end
                        else
                            if xp_fore==xp;
                                xp_fore=xp_fore-1;
                            end
                            d0=sqrt((xp(i)-xp(i+1))^2+(yp(i)-yp(i+1))^2);
                            k=k+P1(xp_fore+d_mov,yp_fore+d_mov)*d0;
                        end
                    end
                    
                end
            end
            
            %%第四种情况，射线过右边
            
            if cross2<0
                if cross4<=0 %与右边相交的情况下与左边相交
                    xp4=zeros(r,1);
                    yp4=zeros(r,1);
                    for i=1:r
                        xp4(i)=r-i+1-d_mov;
                        yp4(i)=tan(gamma)*(xp4(i)-Sr)+Sc;%射线与x=i相交的交点(xp4,yp4)
                    end
                    for i=1:r-1
                        xp4_fore=floor(roundn(xp4(i),-4)+d_mov)-d_mov;
                        yp4_fore=floor(roundn(yp4(i),-4)+d_mov)-d_mov;
                        yp4_back=floor(roundn(yp4(i+1),-4)+d_mov)-d_mov;
                        if yp4_fore==yp4_back %说明相邻两个点在同一个像素区域内
                            d4=sqrt(1+(yp4(i)-yp4(i+1))^2);
                            if xp4_fore==xp;
                                xp4_fore=xp_fore-1;
                            end
                            k=k+P1(xp4_fore+d_mov,yp4_fore+d_mov)*d4;
                        else %说明相邻两个点在相邻两个像素区域内
                            yp4_mid=floor(roundn(yp4(i+1),-4)+d_mov)-d_mov;
                            xp4_mid=(yp4_mid-Sc)/tan(gamma)+Sr;
                            d41=sqrt((xp4(i)-xp4_mid)^2+(yp4(i)- yp4_mid)^2);
                            d42=sqrt((xp4(i+1)-xp4_mid)^2+(yp4(i+1)- yp4_mid)^2);
                            if xp4_fore==xp;
                                xp4_fore=xp_fore-1;
                            end
                            k=k+P1(xp4_fore+d_mov,yp4_fore+d_mov)*d41+P1(xp4_fore+d_mov,yp4_mid+d_mov)*d42;
                        end
                    end
                end
                
                if cross3<0 %与右边相交的情况下与上边相交
                    xp_end=(c_mov-Sc)/tan(gamma)+Sr; %射线与上边交点(xp4,c_mov)
                    xp_end=roundn(xp_end,-4);
                    num=r-floor(xp_end+d_mov)+1; %射线与图像纵轴线相交的点数
                    xp=zeros(num,1);
                    yp=zeros(num,1);
                    %获取与图像纵轴线相交的点所有坐标
                    for i=1:num-1
                        xp(i)=r-i+1-d_mov;
                        yp(i)=tan(gamma)*(xp(i)-Sr)+Sc;
                    end
                    xp(num)=xp_end;
                    yp(num)=c_mov;
                    if num>2
                        for i=1:num-1
                            xp_fore=floor(roundn(xp(i),-4)+d_mov)-d_mov;
                            yp_fore=floor(roundn(yp(i),-4)+d_mov)-d_mov;
                            yp_ceil=ceil(roundn(yp(i+1),-4)+d_mov)-d_mov;
                            dot_num=yp_ceil-yp_fore; %与图像纵轴线相交的相邻两个点之间穿过的像素个数，也等于穿过的射线分段
                            if dot_num>1
                                x_dot=zeros(dot_num-1,1);
                                y_dot=zeros(dot_num-1,1);
                                for j=1:dot_num-1
                                    y_dot(j)=yp_fore+j;%xp(i)与xp(i+1)之间的穿过的像素的点的纵坐标
                                    x_dot(j)=(y_dot(j)-Sc)/tan(gamma)+Sr;%xp(i)与xp(i+1)之间的穿过的像素的点的横坐标
                                end
                                d=zeros(dot_num,1);
                                if dot_num>2
                                    for d_count=2:dot_num-1
                                        d(1)=sqrt((xp(i)-x_dot(1))^2+(yp(i)-y_dot(1))^2);%第i个点与第i+1个点之间穿过像素第一段的距离
                                        d(d_count)=sqrt((x_dot(d_count-1)-x_dot(d_count))^2+(y_dot(d_count-1)-y_dot(d_count))^2);
                                        d(dot_num)=sqrt((x_dot(dot_num-1)-xp(i+1))^2+(y_dot(dot_num-1)-yp(i+1))^2);
                                    end
                                else
                                    d(1)=sqrt((xp(i)-x_dot(1))^2+(yp(i)-y_dot(1))^2);
                                    d(dot_num)=sqrt((x_dot(dot_num-1)-xp(i+1))^2+(y_dot(dot_num-1)-yp(i+1))^2);
                                end
                                
                                for k_count=1:dot_num
                                    if xp_fore==xp;
                                        xp_fore=xp_fore-1;
                                    end
                                    k=k+P1(xp_fore+d_mov,yp_fore+k_count-1+d_mov)*d(k_count);
                                end
                            else
                                if xp_fore==xp;
                                    xp_fore=xp_fore-1;
                                end
                                d0=sqrt((xp(i)-xp(i+1))^2+(yp(i)-yp(i+1))^2);
                                k=k+P1(xp_fore+d_mov,yp_fore+d_mov)*d0;
                            end
                        end
                        
                    end
                end
            end 
        end
        R(s,t+1)=k*200*sqrt(2)*cos(n*p_theta*pi/180);%修正后的投影
        R_VIEW(s,t+1)=k;%投影数据
        
    end
end





for nn=-proj_num/2:proj_num/2
    if(nn==0)
        filterR_L(nn+proj_num/2+1)=1/(8*(p_theta*pi/180)^2);
    elseif(rem(nn,2)==0)
        filterR_L(nn+proj_num/2+1)=0;
    else
        filterR_L(nn+proj_num/2+1)= -1/(2*( pi * sin ( nn * p_theta * pi / 180 ) )^2 ) ;
    end
end
R_temp=zeros(2*(proj_num+1)-1,1);
R_filter=zeros(proj_num+1,c);
for i=1:rot_num+1
    R_temp=conv(R(:,i),filterR_L);
    for j=proj_num+1-proj_num/2:proj_num+1+proj_num/2
        R_filter(j-proj_num+proj_num/2,i)=R_temp(j) * p_theta * pi / 180 ;
    end
end 



%%滤波


%width = 2^nextpow2(size(R,1));  % set width for fft transformation
%proj_fft = fft(R, width);
% Ramp filter function  from 0 to width then to  0




%%反投影
for t=0:rot_num
    Sr=cos(rot_angle*t*pi/180)*Sr0-sin(rot_angle*t*pi/180)*Sc0;
    Sc=sin(rot_angle*t*pi/180)*Sr0+cos(rot_angle*t*pi/180)*Sc0;
    gamma0=45+rot_angle*t+theta/2;%射线左边界
    temp=0;
    rot_anglearch=rot_angle*pi/180;
%     if gamma0>180
%         gamma0=gamma0-180;
%     end
    for x=-r_mov:r_mov
        for y=-c_mov:c_mov
            if (x<=Sr)&&(Sr<=x+1)
                e=(gamma0-90)/p_theta;
                e1=floor(roundn(e,-4));%定位点x，y所处投影线相对位置
                if e>=0&&e<=proj_num
                    temp=R_filter(e1+2,t+1)*(e-e1)+R_filter(e1+1,t+1)*(1-e+e1);%插值
                    temp=rot_anglearch*temp/((x+0.5-Sr)^2+(y+0.5-Sc)^2);%1/L^2加权
                end
            else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%
%                 tanv=(y+0.5-Sc)/(x+0.5-Sr);%点x,y与射线源连线的tan值
                cosv = ( Sr - x - 0.5 ) / sqrt ( ( Sr - x - 0.5 )^2 + ( Sc - y - 0.5 )^2 ) ;
%                 if ( Sc < 0 ) 
%                     Angle = acos ( cosv ) + pi ;
%                 else
%                     Angle = acos ( cosv ) ;
%                 end
                    Angle = acos ( cosv ) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                e=(gamma0-Angle*180/pi)/p_theta;
                e1=floor(roundn(e,-4));%定位点x，y所处投影线相对位置
                if e>=0&&e<=proj_num
                    if roundn(e,-2)~=roundn(e1,-2)
                    temp=R_filter(e1+2,t+1)*(e-e1)+R_filter(e1+1,t+1)*(1-e+e1);%插值
                    else
                        temp=R_filter(e1+1,t+1);%1/L^2加权
                    end
                   temp=rot_anglearch*temp/((x+0.5-Sr)^2+(y+0.5-Sc)^2);%1/L^2加权
                end
            end
            newP1(x+r_mov+1,y+c_mov+1)=newP1(x+r_mov+1,y+c_mov+1)+temp ;
        end
    end
end

% 
% for i = 1:180  
%     % rad is the angle of the projection line , not projection angle  
%     rad = theta(i)*pi/180; 
% 
%     for x = (-128/2+1):128/2  
%         for y = (-128/2+1):128/2  
%             t = round(x*cos(rad+pi/2)+y*sin(rad+pi/2));
% %             if t>-65&&t<65
%             fbp(x+128/2,y+128/2)=fbp(x+128/2,y+128/2)+proj_ifft(t+round(size(R,1)/2),i);  
% %             else
% %                 fbp(x+128/2,y+128/2)=fbp(x+128/2,y+128/2);
%         end  
%     end  
% end  
% fbp = fbp/180;

%p_1=newP1;
%newP1=newP1./(rot_num+1);
%newP1=newP1.*newP1;
%newP1=sqrt(newP1);
figure, imshow(newP1,[]);