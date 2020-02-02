clear all;
close all

width = 7;
high = 2; %2.02 on the picture
w=1; %center conductor width
v=10; %set the center conductor to ten volt
dh=1; %dielectric high
h=0.1; %size of mesh
er = 9.6;
contdh=5;%contour distance, number of meshes away from center conductor
contdv=7;
% matrix size
rmax = high/h-1;
cmax = width/h-1;
tot_node=rmax*cmax; %total number of nodes

Lborder=[1:cmax:tot_node];      %left border
Rborder=[cmax:cmax:tot_node];       %right border
phi_mat=[];         %empty matrix for phi values

for die=1:2
    bc=zeros(tot_node,1);
    
    A=-4*eye(tot_node);
    
    n=1;
    k=1;
    int_line=0;
    mem=zeros(1,2);
  
for node = 1: tot_node
    %find dielectric interface nodes, dh = dielectric high
    if(abs(cell(node/cmax)*h-(high-dh))<eps(1000))
        if(die==2)
            int_line(1,k)=node;
            A(node,node)=-4*(er+1);
            k=k+1;
        end
        %find the center conductor nodes
        if((mod(node,cmax)*h)>=(width-w)/2)
            bc(node,1)=v;
            A(node,node)=1;
            mem(1,n)=node;
            n=n+1;
        end
        
    end
    
    %node on the left of phi0
    if(any(node==Lborder))
    else
        if(die==2 && any(node==int_line))
           A(node,node-1)=1;
        else
            A(node,node-1)=1;
        end
    end
    %node on the top of phi0
    if(node > cmax)
        if(die==2 && any(node==int_line))
            A(node,node-cmax)=2;
        else
            A(node,node-cmax)=1;
        end
        
    end
    %node on the botoom of phi0
    if(node<= (tot_node - cmax))
        if(die==2 && any(node==int_line))
            A(node,node+cmax)=2*er;
        else
            A(node,node+cmax)=1;
        end
    end
    %node on the right of phi0
    if(any(node==Rborder))
        
    else
        if(die==2 && any(node==int_line))
            A(node,node+1)=er+1;
        else
            A(node,node+1)=1;  
 %Lines 84 - 99 Missing
 %following "end" lines are assumed 
        end
    end
end
end
 %{
 
 
 
 
 
 
 
 
 %}
 % %     A((int_line(1,end)+1):end,:)=er*A((int_line(1,end)+1):end,:);
 % %     bc ((int_line(1,end)+1):end,1)=er*bc((int_line(1,end)+1):end,1);
 % % end
 % % take out the known node from the matrix
 % A=[A(:,1:(mem(1,1)-1)),A(:,(mem(1,end)+1):end)];
 % A=[A(1:(mem(1,1)-1),:);A((mem(1,end)+1):end,:)];
 % bc=[bc(1:(mem(1,1)-1,1);bc((mem(1,end)+1):end,1)];
 % % cond =(sum(sum(A.^2)))^0.5*(sum(sum((inv(a)).^2)))^0.5
 for cen=mem(1,1):mem(1,end)
     A(cen,:)=zeros(1,tot_node);
     A(cen,cen)=1;
     
 end
 phi=A\bc;
 % knod=ones(size(mem'))*v;
 % phi=[phi(1:mem(1,1)-1),1);knod;phi(mem(1,1):end,1)];
 % reshape the phi column matix back to geometry dimensions

 for r=1:max
     phi_mat(r,:,die)=phi(1+ (r-1)*cmax : ((r-1)*cmax) +cmax,1);
 end
 %plot phi matrix
 figure;
 imagesc([0:h:width], [0:h:high], phi_mat(:,:,die)); %or pcolor();
 caxis([min(phi(:,:,1)), max(phi(:,:,1))]);colorbar;
 hold on
 %calculate the charge
 pre q=0;
 %set the contour boundry
 hmid=round((high-dh)/h);
 uph=hmid-contdh;
 downh=hmid+condtdh;
 % vmid=cell(width/2/h);
 % leftv= cell(vmid/2);
 % rightv=vmid+leftv;
 
 
