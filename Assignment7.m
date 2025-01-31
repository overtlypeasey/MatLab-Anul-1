N=10;
x=rand(1,N)*10;
y=rand(1,N)*10;
order=[1:N,1];
F=@(x1, y1, x2, y2, x, y) (x-x1)*(y2-y1)-(x2-x1)*(y-y1);
intersect=@(x1, y1, x2, y2, x3, y3, x4, y4) F(x1,y1,x2,y2,x3,y3)*F(x1,y1,x2,y2,x4,y4)<0 && F(x3,y3,x4,y4,x1,y1)*F(x3,y3,x4,y4,x2,y2)<0;
k=0;
no_of_intersection=1;
while no_of_intersection>0

p1=randi(N-1)+1;
p2=randi(N-1)+1;
order([p1 p2])=order([p2 p1]);
    no_of_intersection=0;
for i=1:N-3
    for j=i+2:N
        p1=order(i);
        p2=order(i+1);
        p3=order(j);
        p4=order(j+1);
        if intersect(x(p1),y(p1),x(p2),y(p2),x(p3),y(p3),x(p4),y(p4))
            no_of_intersection=no_of_intersection+1;
           
        end  
    end
end
 

        k=k+1;
end
 plot(x(order),y(order),'-x')
        pause(0.0001)
disp(k);