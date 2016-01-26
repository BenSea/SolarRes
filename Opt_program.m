clear all; close all; clc
format compact;format short g;

% unit conversion prefix
A = 1e-10; nm = 1e-9; um = 1e-6; mm = 1e-3; cm = 1e-2;

% CONSTANTS and other PARAMETERS for the calculation
rho_Au = 2.44e-8; % Resistivity of Gold in [ohm-m] @ 20 degree C
q = 1.602e-19; % [Coulomb]
N_n = (2e17)/(cm^3); 
N_p = (5e18)/(cm^3); 
mu_n = (2.8e3)*(cm^2); % Permeability in [cm^2/(V-sec)]
mu_p = (4e2)*(cm^2); % Permeability in [cm^2/(V-sec)]

% Thickness
t_emit = 500*nm; % Thickness of P-TYPE EMITTER substrate in [nm]
t_base = 500*um; % Thickness of N-TYPE BASE substrate in [um]
t_bm = 2000*A; % Thickness of BACK METAL in [A]
t_fm = 3000*A; % input('front metal thickness (A): ')*A; % [A]

% CELL parameter
L_cell = input('Square cell side (mm): ')*mm; % Cell size in [mm]
A_cell = L_cell^2; % Area of the TPV cell
cvrg = input('Desired coverage area (%): '); % %
A_metal = A_cell * cvrg / 100; % Area of front metal cover
M_cell = 25*um; % input('Cell margin (um): ')*um; % Cell margin in [um]

% BUSBAR parameters
W_bus_m = input('minimum busbar width (um): ')*um; % Minimum busbar width in [um]
W_bus = W_bus_m + (10*um*(0:8*L_cell/mm)); % range of busbar widths 
L_bus = L_cell - (2*M_cell); %- (2*r_pad); % length of the busbar
A_bus = L_bus .* W_bus; % area of the busbar

%FINGER parameters
W_fin_m = 10*um; % Minimum finger width
W_fin = W_fin_m + 2*um*(0:5); % range of finger widths
L_fin = (L_bus - W_bus)/2; % range of finger lengths
A_fin = A_metal - A_bus; % - A_pad % total area of fingers to meet coverage
N_fin = floor((1./W_fin)'*((A_fin/2)./L_fin)); % number of fingers on one side

Vmp = 0.14;
Voc = 0.20;
Isc = 0.0082;
FF = 0.531;

Jmp = (FF*Isc*Voc/Vmp)/(L_cell^2);
rho_ct = 3e-6*(cm^2);


%--CALCULATION
R_fc = 1.66; % front contact resistance 
R_emit = (1/(q*N_p*mu_p)) / t_emit; % 
R_base = (1/(q*N_n*mu_n)) / t_base; % ohm
R_bm = rho_Au / t_bm; % ohm
R_sh = 195;

[M,N] = size(N_fin);
for m = 1:M
    for n = 1:N
        if N_fin(m,n) >= 4
        G_fin(m,n) = um*floor((L_cell-2*M_cell-W_fin(m)*N_fin(m,n))/(N_fin(m,n)-1)/um);
        R_bus = (rho_Au/t_fm)*(G_fin(m,n)/W_bus(n));
        R_fin = (rho_Au/t_fm)*(L_fin(n)/W_fin(m));
        R_L = R_emit*(G_fin(m,n)/L_fin(n));

        %init
        R_A = 0;    R_Y3u = 0;

        R_Dal = R_fin;  R_Dbl = R_bus;  R_Dcl = R_fin + R_L;

        R_Dcru = R_Dcl;

        R_Dbrb = R_fin;
        R_Darb = R_L;

        R_Y1r = R_bus;
        R_Y2r = R_fin;

        for k = 1:N_fin(m,n) - 2;
        %Left side top transformation
            [R_Y1,R_Dcl,R_Y3r] = D2Y(R_Dal,R_Dbl,R_Dcl);
            R_Dbru = R_Y3u + R_Y1;
            R_DcL = R_Dcl + R_L;
        %Right side
            [R_Daru,R_Dbl,R_Dcrb] = Y2D(R_Y1r,R_Y2r,R_Y3r);
            [R_Y1u,R_Y2u,R_Y3u] = D2Y(R_Daru,R_Dbru,R_Dcru);
            R_A = R_A + R_Y1u;
            if k < N_fin(m,n)-2
                [R_Y2r,R_Y2b,R_Y3b] = D2Y(R_Darb,R_Dbrb,R_Dcrb);
                R_Dcru = R_Y2u + R_Y2b;
                R_Darb = R_Y3b + R_L;
            else
                R_sl = R_Dal + R_Dcl;
                R_pl = paraR(R_sl,R_Dbl);
                R_pl = R_pl + R_Y3u;
                R_sr = R_Darb + R_Dbrb;
                R_pr = paraR(R_sr,R_Dcrb);
                R_pr = R_pr + R_Y2u;
                R_p = paraR(R_pl,R_pr);
                R_fm(m,n) = R_A + R_p;
            end
        end        
        else
            G_fin(m,n) = NaN(1);
            R_fm(m,n) = NaN(1);
            N_fin(m,n) = NaN(1);
        end
        P_m(m,n) = (1/12) * L_fin(n)^2 * R_fm(m,n) * Jmp * G_fin(m,n) / (Vmp * W_fin(m)) *100;
    end
end

N_fin = N_fin;
W_bus_um = W_bus / um;

figure()
subplot(2,2,1)
plot3(W_bus_um,N_fin,R_fm,'o-')
xlabel('Busbar Width in \mum')
ylabel('Number of Fingers')
zlabel('Front Metal Resistance in m\Omega')
legend('10','15','20','25','30','35','40','45')
grid on;
subplot(2,2,2);
plot(W_bus_um,R_fm,'*-')
axis([min(W_bus_um) max(W_bus_um) 0 max(max(R_fm))])
xlabel('Busbar Width in \mum')
ylabel('Front Metal Resistance in m\Omega')
legend('10','15','20','25','30','35','40','45')
grid on;
subplot(2,2,3);
plot(N_fin',R_fm','x-')
xlabel('Number of Fingers')
ylabel('Front Metal Resistance in m\Omega')
%legend('10','20','30','40','50','60','70','80')
legend('10','15','20','25','30','35','40','45')
grid on;
subplot(2,2,4);
plot(W_bus_um,N_fin,'.-')
xlabel('Busbar Width in \mum')
ylabel('Number of Fingers')
%legend('10','20','30','40','50','60','70','80')
legend('10','15','20','25','30','35','40','45')
grid on;

figure()
subplot(2,2,1)
plot3(W_bus_um,N_fin,P_m,'o-')
xlabel('Busbar Width in \mum')
ylabel('Number of Fingers')
zlabel('Power Loss due to Metal Layer in %')
legend('10','12','14','16','18','20')
grid on;
subplot(2,2,2);
plot(W_bus_um,P_m,'*-')
axis([min(W_bus_um) max(W_bus_um) 0 max(max(P_m))])
xlabel('Busbar Width in \mum')
ylabel('Power Loss due to Metal Layer in %')
legend('10','15','20','25','30','35','40','45')
grid on;
subplot(2,2,3);
plot(N_fin',P_m','x-')
xlabel('Number of Fingers')
ylabel('Power Loss due to Metal Layer in %')
legend('10','15','20','25','30','35','40','45')
grid on;
subplot(2,2,4);
plot(W_bus_um,N_fin,'.-')
xlabel('Busbar Width in \mum')
ylabel('Number of Fingers')
%legend('10','20','30','40','50','60','70','80')
legend('10','15','20','25','30','35','40','45')
grid on;

[X Y] = size(N_fin);
count = 0;
for y = 1:Y
    for x = 1:X
        count = count + 1;
        result(count,:) = [W_bus_um(y) W_fin(x)/um L_fin(y)/um N_fin(x,y) G_fin(x,y)/um R_fm(x,y)/mm];
    end
end

[x y] = min(R_fm');

for p = 1:length(y)
    result_m(p,:) = [W_bus_um(y(p)) W_fin(p)/um L_fin(y(p))/um N_fin(p,y(p)) G_fin(p,y(p))/um R_fm(p,y(p))/mm];
    result_M(p,:) = [W_bus_um(1) W_fin(p)/um L_fin(1)/um N_fin(p,1) G_fin(p,1)/um R_fm(p,1)/mm];
end

fprintf('\t\tW_Bus[um]\tW_Fin[um]\tL_Fin[um]\t#_of_Fin\tGap_Fin[um]\tR_Series[mohm]\n');
disp(result);

%disp(result_M);
%disp(result_m);

[R,r] = min(result_m(:,6));
%r = 1;
P_emit = R_emit * Jmp * ((result_m(r,5)*um)^2) / (12*Vmp);
P_metal = (1/12) * (result_m(r,3)*um)^2 * (result_m(r,6)*mm) * Jmp * (result_m(r,5)*um) / (Vmp * (result_m(r,2)*um))
P_contact = rho_ct * Jmp * (result_m(r,5)*um) / (Vmp * (result_m(r,2)*um));
Ametal = result_m(r,1)*um*L_bus + result_m(r,2)*result_m(r,3)*um^2*result_m(r,4);
%P_emit = R_emit * Jmp * ((result_M(r,5)*um)^2) / (12*Vmp);
%P_metal = (1/12) * (result_M(r,3)*um)^2 * (result_M(r,6)*mm) * Jmp * (result_M(r,5)*um) / (Vmp * (result_M(r,2)*um));
%P_contact = rho_ct * Jmp * (result_M(r,5)*um) / (Vmp * (result_M(r,2)*um));
%Ametal = result_M(r,1)*um*L_bus + result_M(r,2)*result_M(r,3)*um^2*result_M(r,4);
P_shadow = Ametal/A_cell

P_loss = P_emit + P_metal + P_contact + P_shadow

OptB_W = (L_bus^2)*sqrt((rho_Au/t_fm)*Jmp/(3*Vmp))/um
