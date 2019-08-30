%Brownian Motion with Periodic boundary conditions

N = 25;                             %Number of particles in system
R = 0.2;                            %Radii of the particles
L = 3.0;                            %Size of periodic domain
max_vel = 15.0;                     %Max velocity components of the particles
max_t = 40;                         %Max simulation time
dt = 1e-2;                          %Time step size
m = 1.0;                            %Particle mass

writevideo = true;                  %Generate video of simulation
show_mirror_domains = true;         %Show mirror domains

collisions_in_main_domain_counter = 0;
collisions_with_mirror_domain_counter = 0;

%Initialization
E_wall = L;
W_wall = -L;
N_wall = L;
S_wall = -L;

E_lab = 1;
W_lab = 2;
N_lab = 3;
S_lab = 4;

NE_lab = 5;
NW_lab = 6;
SW_lab = 7;
SE_lab = 8;

C_lab = 9;

err = 1e-8;                         %Small number
h = L - R - err;                    %Domain limit particle injection

X = zeros(N,9);
Y = zeros(N,9);
    
V_X = zeros(N,9);
V_Y = zeros(N,9);

x = zeros(N,1);
y = zeros(N,1);

avg_kin_energy = 0.0;

%Create video
frameps = 32;                       %Set framerate (fps) for video
if writevideo == true
    writerObj = VideoWriter('C:\Users\d-w-h\Desktop\Home\Sim_vid_EC_PBC_with_mirror_long_store_2.avi','Motion JPEG AVI');
    writerObj.FrameRate = frameps;
    open(writerObj);
end

%Generate initial velocities
v_x = 2*max_vel*rand(N,1) - max_vel;
v_y = 2*max_vel*rand(N,1) - max_vel;

%Generating initial positions of particles
for n = 1:N
    x(n) = 2 * h * rand - h;
    y(n) = 2 * h * rand - h;
    overlap = true;
    while overlap == true
        overlap = false;
        for i = 1:(n-1)
            dx = x(i) - x(n);
            dy = y(i) - y(n);
            if dx*dx + dy*dy < 1.0001*(2*R)*(2*R)
                overlap = true;
            end
        end
        if overlap == true
            x(n) = 2 * h * rand - h;
            y(n) = 2 * h * rand - h;
        end        
    end
end

%Check kinetic energy before simulation
%Check average kinetic energy
avg_kin_energy_pre = 0.0;
for n = 1:N
    avg_kin_energy_pre = avg_kin_energy_pre + 0.5 * m * (v_x(n) * v_x(n) + v_y(n) * v_y(n));
end
avg_kin_energy_pre = avg_kin_energy_pre / N;

%Start simulation code
time = 0;
frame_counter = 0;
while time < max_t   
    %Copying velocities to mirror domains
    for lab = E_lab:C_lab
        V_X(:,lab) = v_x;
        V_Y(:,lab) = v_y;
    end

    %Copying locations to mirror domains
    X(:,C_lab) = x;
    Y(:,C_lab) = y;

    X(:,E_lab) = x + 2 * L;
    Y(:,E_lab) = y;
    X(:,W_lab) = x - 2 * L;
    Y(:,W_lab) = y;
    X(:,N_lab) = x;
    Y(:,N_lab) = y + 2 * L;
    X(:,S_lab) = x;
    Y(:,S_lab) = y - 2 * L;

    X(:,NE_lab) = x + 2 * L;
    Y(:,NE_lab) = y + 2 * L;
    X(:,NW_lab) = x - 2 * L;
    Y(:,NW_lab) = y + 2 * L;
    X(:,SE_lab) = x + 2 * L;
    Y(:,SE_lab) = y - 2 * L;
    X(:,SW_lab) = x - 2 * L;
    Y(:,SW_lab) = y - 2 * L;

    %Initialise collision related info
    coll_partner_1 = 0;
    coll_partner_2 = 0;
    coll_outside_domain = false;
    domain_lab = 0;
    coll_time = 1e+10;

    %Checking collision time between particles in main domain and mirror domains
    for n = 1:N
        for lab = E_lab:SE_lab
            for i = 1:N
                if i ~= n
                    rab = [x(n) - X(i, lab);y(n) - Y(i, lab)];
                    vab = [v_x(n) - V_X(i, lab);v_y(n) - V_Y(i, lab)];
                    Disc = (rab'*vab)^2-vab'*vab*(rab'*rab-(2*R)^2);
                    coll_time_particle = dt + err;

                    if Disc > 0
                        coll_time_particle = ((-rab'*vab) - sqrt(Disc)) / (vab'*vab);
                    end

                    if coll_time_particle < coll_time && coll_time_particle >= 0
                        coll_time = coll_time_particle;
                        coll_outside_domain = true;
                        domain_lab = lab;
                        coll_partner_1 = n;
                        coll_partner_2 = i;

                    end
                end
            end
        end
    end

    %Checking collision time between particles within the main domain
    for n = 1:N
        for i = (n+1):N
            rab = [x(n) - x(i);y(n) - y(i)];
            vab = [v_x(n) - v_x(i);v_y(n) - v_y(i)];
            Disc = (rab'*vab)^2-vab'*vab*(rab'*rab-(2*R)^2);
            coll_time_particle = dt + err;

            if Disc > 0
                coll_time_particle = ((-rab'*vab) - sqrt(Disc)) / (vab'*vab);
            end

            if coll_time_particle < coll_time && coll_time_particle >= 0
                coll_time = coll_time_particle;
                coll_outside_domain = false;
                domain_lab = 0;
                coll_partner_1 = n;
                coll_partner_2 = i;
            end
        end
    end

    %Update positions and velocities 
    if coll_time < dt 
        %Update positions to the point of collision minus some small, negligible,  numerical value to
        %prevent particles from getting stuck in a wall
        x = x + v_x*coll_time*(1-err);
        y = y + v_y*coll_time*(1-err);
        
        %Handle periodic boundaries
        for n = 1:N
            if x(n) < W_wall
                x(n) = x(n) + 2 * L;
            elseif x(n) > E_wall
                x(n) = x(n) - 2 * L;
            end
            if y(n) < S_wall
                y(n) = y(n) + 2 * L;
            elseif y(n) > N_wall
                y(n) = y(n) - 2 * L;
            end
        end

        time = time + coll_time;
           
        %Update mirror particle positions
        X(:,C_lab) = x;
        Y(:,C_lab) = y;

        X(:,E_lab) = x + 2 * L;
        Y(:,E_lab) = y;
        X(:,W_lab) = x - 2 * L;
        Y(:,W_lab) = y;
        X(:,N_lab) = x;
        Y(:,N_lab) = y + 2 * L;
        X(:,S_lab) = x;
        Y(:,S_lab) = y - 2 * L;

        X(:,NE_lab) = x + 2 * L;
        Y(:,NE_lab) = y + 2 * L;
        X(:,NW_lab) = x - 2 * L;
        Y(:,NW_lab) = y + 2 * L;
        X(:,SE_lab) = x + 2 * L;
        Y(:,SE_lab) = y - 2 * L;
        X(:,SW_lab) = x - 2 * L;
        Y(:,SW_lab) = y - 2 * L;
                
        if coll_outside_domain == false
            collisions_in_main_domain_counter = collisions_in_main_domain_counter + 1;
            %Update velocities of colliding particles
            ra = [x(coll_partner_1);y(coll_partner_1)];
            va = [v_x(coll_partner_1);v_y(coll_partner_1)];
            rb = [x(coll_partner_2);y(coll_partner_2)];
            vb = [v_x(coll_partner_2);v_y(coll_partner_2)];              
            n = (ra-rb)/sqrt((ra-rb)'*(ra-rb));

            del_v = n*((va - vb)'*n);

            v_x(coll_partner_1) = v_x(coll_partner_1) - del_v(1);
            v_y(coll_partner_1) = v_y(coll_partner_1) - del_v(2);
            v_x(coll_partner_2) = v_x(coll_partner_2) + del_v(1);
            v_y(coll_partner_2) = v_y(coll_partner_2) + del_v(2);            
            
        elseif coll_outside_domain == true
            collisions_with_mirror_domain_counter = collisions_with_mirror_domain_counter + 1;           
            %Update velocities of colliding particles
            ra = [x(coll_partner_1);y(coll_partner_1)];
            va = [v_x(coll_partner_1);v_y(coll_partner_1)];
            rb = [X(coll_partner_2, domain_lab);Y(coll_partner_2, domain_lab)];
            vb = [V_X(coll_partner_2, domain_lab);V_Y(coll_partner_2, domain_lab)];              
            n = (ra-rb)/sqrt((ra-rb)'*(ra-rb));
            
            del_v = n*((va - vb)'*n);
            
            %Instead of updating V_X of particle n in the mirror domain
            %and then updating the velocity in the real domain the
            %velocities of n in the real domain are computed directly                      
            v_x(coll_partner_1) = v_x(coll_partner_1) - del_v(1);
            v_y(coll_partner_1) = v_y(coll_partner_1) - del_v(2);
            v_x(coll_partner_2) = v_x(coll_partner_2) + del_v(1);
            v_y(coll_partner_2) = v_y(coll_partner_2) + del_v(2); 
        end
    elseif coll_time >= dt %Update positions and velocities using dt
        x = x + v_x*dt;
        y = y + v_y*dt;

        %Handle periodic boundaries
        for n = 1:N
            if x(n) < W_wall
                x(n) = x(n) + 2 * L;
            elseif x(n) > E_wall
                x(n) = x(n) - 2 * L;
            end
            if y(n) < S_wall
                y(n) = y(n) + 2 * L;
            elseif y(n) > N_wall
                y(n) = y(n) - 2 * L;
            end
        end
        
        %Update mirror positions for figure
        X(:,E_lab) = x + 2 * L;
        Y(:,E_lab) = y;
        X(:,W_lab) = x - 2 * L;
        Y(:,W_lab) = y;
        X(:,N_lab) = x;
        Y(:,N_lab) = y + 2 * L;
        X(:,S_lab) = x;
        Y(:,S_lab) = y - 2 * L;

        X(:,NE_lab) = x + 2 * L;
        Y(:,NE_lab) = y + 2 * L;
        X(:,NW_lab) = x - 2 * L;
        Y(:,NW_lab) = y + 2 * L;
        X(:,SE_lab) = x + 2 * L;
        Y(:,SE_lab) = y - 2 * L;
        X(:,SW_lab) = x - 2 * L;
        Y(:,SW_lab) = y - 2 * L;        
        
        time = time + dt;
    end
    
    TIME = time
    
    if frame_counter == floor(time/dt)
        frame_counter = frame_counter + 1;

        if show_mirror_domains == true
            scatter(x,y, 2 * R * 200, 'b', 'filled')
            
            hold on
            scatter(X(:,E_lab),Y(:,E_lab), 2 * R * 200, 'k', 'filled')
            scatter(X(:,W_lab),Y(:,W_lab), 2 * R * 200, 'k', 'filled')
            scatter(X(:,N_lab),Y(:,N_lab), 2 * R * 200, 'k', 'filled')
            scatter(X(:,S_lab),Y(:,S_lab), 2 * R * 200, 'k', 'filled')

            scatter(X(:,NE_lab),Y(:,NE_lab), 2 * R * 200, 'k', 'filled')
            scatter(X(:,NW_lab),Y(:,NW_lab), 2 * R * 200, 'k', 'filled')
            scatter(X(:,SE_lab),Y(:,SE_lab), 2 * R * 200, 'k', 'filled')
            scatter(X(:,SW_lab),Y(:,SW_lab), 2 * R * 200, 'k', 'filled')
            
            xlim([(-L - 2 * L) (L + 2 * L)])
            ylim([(-L - 2 * L) (L + 2 * L)])            
            hold off
        else
            scatter(x,y, 2 * R * 800, 'b', 'filled')
            xlim([(-L) (L)])
            ylim([(-L) (L)])              
        end 
        
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer')  
        if writevideo == 1
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
        end
    end
    
    %Check average kinetic energy
    avg_kin_energy = 0.0;
    for n = 1:N
        avg_kin_energy = avg_kin_energy + 0.5 * m * (v_x(n) * v_x(n) + v_y(n) * v_y(n));
    end
    avg_kin_energy = avg_kin_energy / N
    
    %Check if particles are within bounds
    in_boundaries = false;
    for n = 1:N
        in_boundaries = (x(n) <= E_wall) && (x(n) >= W_wall && ... 
                         y(n) <= N_wall) && (y(n) >= S_wall);
        if in_boundaries
            in_boundaries = true;
        end
    end
    in_boundaries;
    
    %Check overlap
    there_is_overlap = false;
    for n = 1:N
        for i = (n+1):N 
            dx = x(i) - x(n);
            dy = y(i) - y(n);
            if dx*dx + dy*dy < (2*R)*(2*R)
                there_is_overlap = true;
            end    
        end
    end
    there_is_overlap;
end

collisions_in_main_domain_counter
collisions_with_mirror_domain_counter
    
KIN_ENERGY_RATIO = avg_kin_energy_pre / avg_kin_energy
avg_kin_energy_pre

if writevideo == true
    close(writerObj);
end