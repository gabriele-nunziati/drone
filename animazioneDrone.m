%%GRUPPO S PROGETTO 1A
%%Stefano Favale, Gabriele Nunziati, Giorgia Quarta, Jacopo Jop
%% animation script
T=0.01;
Y=step(W*LL/(1+LL),5);
fhand=figure(20);
clf;
max_T=max(T);
min_Y=min(Y);
max_Y=max(Y);
Y_ss=Y(end);
axis([-10 10  -10 10])



%vertici
AA=[-5, -1];
BB=[5, -1];
CC=[5, 1];
DD=[-5, 1];

for i=1:1:500
figure(20);
clf;
    axis([-10 10  -10 10])
    % Plot one frame...
    
    theta= -(pi/4 + Y(i));%angolo
    hold on;
    
    r=[cos(theta) -sin(theta); sin(theta) cos(theta)];%matrice rotazione
    
    
    A1=AA*r;
    B1=BB*r;
    C1=CC*r;
    D1=DD*r;
    
    
    patch( [A1(1) B1(1) C1(1) D1(1)] ,[A1(2) B1(2) C1(2) D1(2)],'b') % 
    
    
    
    axis([-10 10  -10 10])
    title('Risposta al gradino W 1(t)')
    if i==1
        pause(1);
    end
    
    pause(0.01);
    if ishandle(fhand)==false; break;
    end
end