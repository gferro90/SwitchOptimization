clc 
clear all
close all
%% DEFINIZIONE TEMPO DI SIMULAZIONE

% ottimizzazione offline

% Definizione di un predefinito numero di campioni spalmati sull'intervallo
% di tempo
numeroCampioni=20000;

% Definizione del tempo finale di simulazione
tempoFinale=2;

% Definizione del delta t
delta_t=tempoFinale/numeroCampioni;

numeroCampioni=numeroCampioni+1;

% Definizione del vettore dei tempi
t=[0:delta_t:tempoFinale];

%% INPUT DATI
% Input da questo file.
code=0;

if(code==1)
    %inserire il numero di salti:
    numControllori=input('Inserire il numero desiderato di Controllori: ');
    numSalti=numControllori-1;
    
    %inserire il riferimento (quello a regime permanente):
    numRifStatico=input('Inserire il numeratore del riferimento statico: ');
    denRifStatico=input('Inserire il denominatore del riferimento statico: ');
    rifStaticoTempo=input('Inserire l espressione del riferimento statico nel tempo: ');
    
    %inserire il riferimento dinamico:
    numRifDinamico=input('Inserire il numeratore del riferimento dinamico: ');
    denRifDinamico=input('Inserire il denominatore del riferimento dinamico: ');
    rifDinamicoTempo=input('Inserire l espressione del riferimento dinamico nel tempo: ');

    
    %inserire il processo:
    numProcesso=input('Inserire il numeratore del processo: ');
    
    %inserire condizioni iniziali del processo:
    condInizialiProcesso=input('Inserire le condizioni iniziali del processo: ');
    
    for i=1:numControllori
        stringaControllerNum=strcat('Inserire il numeratore del controllore ',num2str(i),': ');
        controllerNum=input(stringaControllerNum);
        stringaControllerDen=strcat('Inserire il denominatore del controllore ',num2str(i),': ');
        controllerDen=input(stringaControllerDen);
        sistema(i).fdtControllore=tf(controllerNum, controllerDen);
    end
    
    guadagnoErrore=input('Inserire guadagno relativo all errore: ');
    guadagnoDerivataErrore=input('Inserire guadagno relativo alla derivata dell errore');
    
else
    numSalti=1;
    numControllori=numSalti+1;

    esponente=500;
    
     numRifDinamico=esponente;
     denRifDinamico=[1 esponente 0];
%     numRifDinamico=1;
%     denRifDinamico=[1 0];
    rifDinamicoTempo=1-exp(-esponente*t);
    
    wn=2*pi/0.3;
    delta=0.8;
    % P=k/(s^2+as+b)
    k=36*(wn^2);
    a=2*delta*wn;
    b=wn^2;
    
    numProcesso=[k];
    denProcesso=[1 a b];
    
    CI_Processo=[0; 0];
    CI_Controllore=0;
    CI_cc=[CI_Controllore; CI_Processo];
    
    Kp=0.1;
    %Ki=[0.2 3];
   Ki=[0.2 3 0.2];

    sistema(1).fdtControllore=tf([Kp Ki(1)], [1 0]);
    sistema(2).fdtControllore=tf([Kp Ki(2)], [1 0]);
    sistema(3).fdtControllore=tf([Kp Ki(3)], [1 0]);

    % guadagno associato all'errore nel funzionale di costo
    p1=10;
    % guadagno associato alla derivata dell'errore nel funzionale di costo
    p2=1e-2;
end

%% COSTRUZIONE FUNZIONI DI TRASFERIMENTO

%fdt del riferimento dinamico
fdtRifDinamico=tf(numRifDinamico, denRifDinamico)

%fdt del processo
fdtProcesso=tf(numProcesso, denProcesso)

%fdt dei sistemi a ciclo chiuso 
for i=1:numControllori
    %il sistema a ciclo chiuso
    sistema(i).fdtCicloChiuso=feedback(sistema(i).fdtControllore*fdtProcesso,1)
end


%% COSTRUZIONE MATRICI DI STATO
% phi = rd - rs
% COSTRUZIONE A_phi
   
% Passare al sistema in forma minima per semplificare la funzione di
% trasferimento
sistema_rd=ss(fdtRifDinamico,'minimal');
   
% Ottenere numeratore e denominatore della fdt
[num_rd,den_rd]=tfdata(tf(sistema_rd),'v');
    
% Calcolo del numero di stati
numStati_r=length(den_rd)-1;

A_rd=zeros(numStati_r, numStati_r);

% Calcolo della matrice A del sistema
if(numStati_r>1)
    A_rd(1:numStati_r-1,2:numStati_r)=eye(numStati_r-1);
    A_rd(1:numStati_r-1,1)=zeros(numStati_r-1,1);
    A_rd(numStati_r,:)=-den_rd(numStati_r+1:-1:2);
else
    if (numStati_r==1)
        A_rd=-den_rd(2);
    end
end
 
% Costruzione del sistema processo
[num_p,den_p]=tfdata(fdtProcesso,'v');
numStati_p=length(den_p)-1;
if(numStati_p>1)
    % 1) Riempire le prime n-1 righe dalla seconda colonna in poi con I shiftato a destra di 1.
    A_p(1:numStati_p-1,2:numStati_p)=eye(numStati_p-1);
    % 2) Riempire la prima colonna fino a n-1 righe di zeri.
    A_p(1:numStati_p-1,1)=zeros(numStati_p-1,1);
    % 3) Riempire l'ultima riga con i coefficienti del denominatore
    % cambiati di segno.
    A_p(numStati_p,:)=-den_p(numStati_p+1:-1:2);
else
    % In questo caso A è scalare.
    A_p=-den_p(2);
end

% Costruire B in forma canonica di controllore:
% 1) Riempire le prime n-1 righe di zeri.
B_p=zeros(numStati_p,1);
% 2) Settare a 1 l'ultima riga.
B_p(numStati_p,1)=1;

% Costruire C in forma canonica di controllore:
C_p=zeros(1,numStati_p);
% 1) Riempire la prima riga di C con i coefficienti del numeratore dal
% grado + basso al + alto.
C_p(1:numStati_p)=num_p(numStati_p+ 1:-1:2);

% Normalizzo la matrice C mettendo i suoi coefficienti nella dinamica
mul=eye(numStati_p);
for j=1:numStati_p
    if(C_p(j)~=0)
        mul(j,j)=C_p(j);
        C_p(j)=1;
    end
end
A_p=mul*A_p*(mul^-1);
B_p=mul*B_p;

Ifull=eye(numStati_r);
Ir=Ifull(1,:);



for i=1:numControllori
     
    % Ricava numeratore e denominatore della fdt a ciclo chiuso. 
    % I vettori sono di uguale lunghezza e vanno dal grado maggiore a
    % quello minore.
    [num_cc, den_cc]=tfdata(sistema(i).fdtCicloChiuso,'v')
    
    % Ricava il numero di stati a ciclo chiuso uguale al grado del
    % denominatore della fdt.
    numStati_cc=length(den_cc)-1;
    
    % Costruire Acc considerando come stato [xc xp]
    % Costruzione del sistema controllore 
    [num_c,den_c]=tfdata(sistema(i).fdtControllore,'v');

    numStati_c=length(den_c)-1;
    A_c=zeros(numStati_c, numStati_c);
    if(numStati_c>1)
        % 1) Riempire le prime n-1 righe dalla seconda colonna in poi con I shiftato a destra di 1.
        A_c(1:numStati_c-1,2:numStati_c)=eye(numStati_c-1);
        % 2) Riempire la prima colonna fino a n-1 righe di zeri.
        A_c(1:numStati_c-1,1)=zeros(numStati_c-1,1);
        % 3) Riempire l'ultima riga con i coefficienti del denominatore
        % cambiati di segno.
        A_c(numStati_c,:)=-den_c(numStati_c+1:-1:2);
    else
        % In questo caso A è scalare.
        A_c=-den_c(2);
    end

    
    % Costruire B in forma canonica di controllore:
    % 1) Riempire le prime n-1 righe di zeri.
    B_c=zeros(numStati_c,1);
    % 2) Settare a 1 l'ultima riga.
    B_c(numStati_c,1)=1;
    
    % Costruire C in forma canonica di controllore:
    C_c=zeros(1,numStati_c);
    % 1) Riempire la prima riga di C con i coefficienti del numeratore dal
    % grado + basso al + alto.
    C_c(1:numStati_c)=num_c(numStati_c+ 1:-1:2);
    
    % Se != 0 D è il coefficiente del termine di grado massimo
    D_c=num_c(1);
           
    % Normalizzo la matrice C mettendo i suoi coefficienti nella dinamica
    mul=eye(numStati_c);
    for j=1:numStati_c
        if(C_c(j)~=0)
            mul(j,j)=C_c(j);
            C_c(j)=1;
        end
    end
    A_c=mul*A_c*(mul^-1);
    B_c=mul*B_c;
    
    
    sistema(i).A_c=A_c;
    sistema(i).B_c=B_c;
    sistema(i).C_c=C_c;
    sistema(i).D_c=D_c;
   
    
    % Costruzione sistema a ciclo chiuso
    sistema(i).A_cc=[A_c -B_c*C_p; B_p*C_c A_p-B_p*D_c*C_p];
    sistema(i).B_cc=[B_c; B_p*D_c];
    sistema(i).C_cc=[zeros(1, numStati_c) C_p];
    sistema(i).D_cc=0;
    
    
    
    if(i==1)
        % La risposta libera del sistema a ciclo chiuso nel dominio di Laplace
        fdtRispostaLiberaCicloChiuso=tf(ss(sistema(i).A_cc, CI_cc, sistema(i).C_cc, sistema(i).D_cc));
    else
        fdtRispostaLiberaCicloChiuso=0;
    end
    
    % L'errore in evoluzione libera
    sistema(i).fdt_ed=(1-sistema(i).fdtCicloChiuso)*fdtRifDinamico-fdtRispostaLiberaCicloChiuso;

    % Passare al sistema in forma minima per semplificare la funzione di
    % trasferimento
    sistema_ed=ss(sistema(i).fdt_ed,'minimal');

    % Ottenere numeratore e denominatore della funzione in Laplace di y-rs
    [num_ed, den_ed]=tfdata(tf(sistema_ed), 'v');
    
    % Calcolare la dimensione della matrice Ai
    numStati_ed=length(den_ed)-1;
    
    % Salvare i dati del primo sistema che serviranno poi per il calcolo
    % delle condizioni iniziali
    if(i==1)
        num_ed_i=num_ed;
        den_ed_i=den_ed;
        numStati_ed_i=numStati_ed;
    end
    
    % Costruzione della matrice Ai del sistema in evoluzione libera
    % ded=Ai*ed con ed=rd-y
    sistema(i).A_ed=zeros(numStati_ed, numStati_ed);
    
    if(numStati_ed>1)
        sistema(i).A_ed(1:numStati_ed-1, 2:numStati_ed)=eye(numStati_ed-1);
        sistema(i).A_ed(1:numStati_ed-1, 1)=zeros(numStati_ed-1,1);
        sistema(i).A_ed(numStati_ed,:)=-den_ed(numStati_ed+1:-1:2);
    else
        sistema(i).A_ed=-den_ed(2);
    end

    numStati=numStati_ed;
    
    % calcolo matrice x il successivo cambio di coordinate a [e e' e'' ...]
    sistema(i).M=zeros(numStati, numStati);
    sistema(i).O=zeros(numStati, numStati_cc);
    
    % chain_1 is A^n    
    chain_1=eye(numStati_cc);
    % chain_2 is [A^(n-1)B | A^(n-2)B |...| B]
    chain_2=zeros(numStati_cc,numStati_cc);
    for j=1:numStati
        sistema(i).O(j,1:numStati_cc)=sistema(i).C_cc*(chain_1);
        sistema(i).M(j,1:numStati_cc)=sistema(i).M(j,1:numStati_cc)+sistema(i).C_cc*(chain_2);
        if(j<numStati)
            chain_1=sistema(i).A_cc*chain_1;
            chain_2=sistema(i).A_cc*chain_2;
            chain_2(:,j)=sistema(i).B_cc;
        end        
    end
    %sistema(i).M=eye(numStati, numStati)-sistema(i).M;
    % Mettere I per lo stato dell' "errore" dinamico ex_rd = rs-rd
    % alla fine facendo [x_r | 0] - M [x_cc | x_r | ex_rd] = [e de ... ex_rd]
    
       
    % numero stati iniziali
    numStati_i=numStati_ed_i;

    % Definizione del parametro epsilon per pesare gli stati non necessari
    epsilon=0;
    
    % Definizione della matrice Q usata nel funzionale di costo
    %sistema(i).Q=zeros(numStati, numStati);
    sistema(i).Q=eye(numStati, numStati)*epsilon;

    
    % Inserire i guadagni che pesano l'errore dinamico e la sua derivata
    sistema(i).Q(1,1)=p1;
    sistema(i).Q(2,2)=p2;
    
    % Calcolo della matrice P utilizzata nella legge di controllo switch
    sistema(i).P=lyap((sistema(i).A_ed)',sistema(i).Q);
        
    
    N=zeros(numStati, numStati_r);
    for kk=1:numStati
       N(kk,:)=Ir*(A_rd^(kk-1));
    end
    
    % Calcolo della matrice T per questo sistema
    sistema(i).T=[-sistema(i).O (eye(numStati)-sistema(i).M)*N];
    sistema(i).Tc=sistema(i).T(:,1:numStati_c);
    sistema(i).Tp=sistema(i).T(:,numStati_c+1:numStati_c+numStati_p);
    sistema(i).Tr=sistema(i).T(:,numStati_c+numStati_p+1:numStati_c+numStati_p+numStati_r);
    
    if(i>1)
        sistema(i-1).Kc=sistema(i-1).Tc'*sistema(i).P*sistema(i-1).Tc;
        sistema(i-1).Kp=sistema(i-1).Tc'*sistema(i).P*sistema(i-1).Tp;
        sistema(i-1).Kr=sistema(i-1).Tc'*sistema(i).P*sistema(i-1).Tr;
    end
    
%     if(i>1)
%         % Calcolo matrici per cambio di coordinate da xcc a xe
%         sistema(i-1).F=sistema(i-1).O'*sistema(i).P*sistema(i-1).O;
%         sistema(i-1).G=sistema(i-1).O'*sistema(i).P*sistema(i-1).M;
%     end
end



%% CALCOLO DELLE CONDIZIONI INIZIALI DEL SISTEMA

% Per Calcolare le condizioni iniziali del sistema ci avvaliamo di un
% semplice metodo per il quale M.ci=R. L'equazione è:
% [b1 b2 .... b_n] [f(0)      ]    [  a0  ]
% [b2 b3 .... 0  ] [f'(0)     ] =  [  a1  ]
% [...        ...] [...       ]    [ ...  ]
% [b_n   .....0  ] [df_n-1(0) ]    [a_n-1 ]
% con b_i coefficienti del denominatore e a_i coefficienti del numeratore.

% CALCOLO CONDIZIONI INIZIALI RIFERIMENTO DINAMICO

% Calcolo della matrice M
M_ParametriCI_rd=zeros(numStati_r);

% Costruzione del vettore dei termini noti [a0 a1 ... a_n-1]
R_TerminiNotiCI_rd=(num_rd(numStati_r+1:-1:2))';

% Construzione della matrice dei parametri
for i=1:numStati_r
    % i coefficienti del denominatore vanno da numStati a 1
    M_ParametriCI_rd(i,1:numStati_r+1-i)=den_rd(numStati_r+1-i:-1:1)';
end
CI_rd=zeros(numStati_r,1);
% Calcolo delle condizioni iniziali come vettore colonna
CI_rd(1:numStati_r,1)=(M_ParametriCI_rd^-1)*R_TerminiNotiCI_rd
% for kk=numStati_r+1:numStati
%     CI_rd(kk,1)=N(kk-1, :)*CI_rd;
% end


% CALCOLO CONDIZIONI INIZIALI DEL SISTEMA ed=rd-y

% Calcolo della matrice M
M_ParametriCI_ed=zeros(numStati_ed_i);

% Costruzione del vettore dei termini noti [a0 a1 ... a_n-1]
R_TerminiNotiCI_ed=(num_ed_i(numStati_ed_i+1:-1:2))';

% Construzione della matrice dei parametri
for i=1:numStati_ed_i
    % i coefficienti del denominatore vanno da numStati a 1
    M_ParametriCI_ed(i,1:numStati_ed_i+1-i)=den_ed_i(numStati_ed_i+1-i:-1:1)';
end

% Calcolo delle condizioni iniziali come vettore colonna
CI_ed=(M_ParametriCI_ed^-1)*R_TerminiNotiCI_ed


%% COSTRUZIONE MATRICI PER OTTIMIZZAZIONE

for kk=1:numSalti
    Mc(numStati*(kk-1)+1:numStati*kk,numStati*(kk-1)+1:numStati*kk)=sistema(kk+1).P-sistema(kk).P; 
end



timeMesh = zeros(numSalti, numeroCampioni);
Aev=zeros(numStati, numStati*numSalti);
for kk=1:numSalti
    Aev(:, numStati*(kk-1)+1:kk*numStati)=sistema(kk).A_ed;
    time=0;
    for kkk=1:numeroCampioni
        timeMesh(kk, kkk)=time;
        time=time+delta_t;
    end
end

cnt=1;
Jnow=Inf;
indexes=zeros(1, numSalti);



%% PLOT RISULTATI
if(numSalti==1)
    Jstore=zeros(1, numeroCampioni);
    [index_opt, Jmin, Jstore]=FindJCostminStore(cnt, indexes, numSalti,numStati, numeroCampioni, timeMesh, Mc, Aev, CI_ed, Jnow, indexes,Jstore);
    J0=ones(1,length(Jstore))*(CI_ed'*sistema(1).P*CI_ed);
    figure(1);
    plot(t,Jstore+J0);
    
    Jmin
    t(index_opt(1))
else
    if(numSalti==2)
        Jstore=zeros(1, numeroCampioni^2);
        [index_opt, Jmin, Jstore]=FindJCostminStore(cnt, indexes, numSalti,numStati, numeroCampioni, timeMesh, Mc, Aev, CI_ed, Jnow, indexes, Jstore);
        Jmatrix=zeros(numeroCampioni, numeroCampioni);
        for i=1:numeroCampioni
            for j=1:numeroCampioni
                Jmatrix(i,j)=Jstore((i-1)*numeroCampioni+j);
            end
        end
        mesh(t,t,Jmatrix);
        
    else
        
        [index_opt, Jmin]=FindJCostmin(cnt, indexes, numSalti,numStati, numeroCampioni, timeMesh, Mc, Aev, CI_ed, Jnow, indexes)
        
        istante_switch=t(index_opt)
    end
end
samples=200000;
% Definizione del delta t
dt=tempoFinale/samples;

samples=samples+1;

% Definizione del vettore dei tempi
t1=[0:dt:tempoFinale];

% Condizione iniziale dello stato del riferimento
x_rd=zeros(numStati_r,1);
x_rd=CI_rd;

% evoluzione dei sistemi autonomi
d_x_cc_s=zeros(numStati_cc,numControllori);
x_cc_s=zeros(numStati_cc,numControllori);
for nSist=1:numControllori
    x_cc_s(:,nSist)=CI_cc;
end
x_cc_system_plot=zeros(numControllori, samples);
x_rd_plot=zeros(numStati_r, samples);
x_ed=zeros(numStati, 1);
x_ed_plot=zeros(numStati, samples);
x_ed(:,1)=CI_ed;

x_cc_normal=CI_cc;
x_cc_normal_plot=zeros(1,samples);
index_normal=1;
switch_i=1;



for i=1:samples
    x_cc_normal_plot(i)=sistema(index_normal).C_cc*x_cc_normal;
    x_rd_plot(:,i)=x_rd;
    x_ed_plot(:,i)=x_ed;
    
    % Evoluzione dei sistemi 
    % evolutione dello stato senza ottimizzazione con switch al tempo
    % ottimo
    if(switch_i<=numSalti)
        if(t1(i)>=t(index_opt(switch_i)))
            index_normal=index_normal+1;
            
            switch_i=switch_i+1;
        end
    end
    
    d_x_cc_normal=sistema(index_normal).A_cc*x_cc_normal+sistema(index_normal).B_cc*x_rd(1,1);
    x_cc_normal=x_cc_normal+d_x_cc_normal*dt;
    
    for nSist=1:numControllori
        d_x_cc_s(:,nSist)=sistema(nSist).A_cc* x_cc_s(:,nSist)+sistema(nSist).B_cc*x_rd(1,1);
        x_cc_s(:,nSist)=x_cc_s(:,nSist)+d_x_cc_s(:,nSist)*dt;
        x_cc_system_plot(nSist,i)=sistema(nSist).C_cc*x_cc_s(:,nSist);
    end
    % Evoluzione di rd
    d_x_rd=A_rd*x_rd;
    x_rd=x_rd+d_x_rd*dt;
    
    d_x_ed=sistema(index_normal).A_ed*x_ed;
    x_ed=x_ed+d_x_ed*dt;
end

Plot1=figure(2);
%ax(1)=subplot(2,1,1);
plot(t1,x_cc_normal_plot,t1,x_rd_plot(1,:), t1, x_cc_system_plot(1,:), t1, x_cc_system_plot(2,:), t1, x_ed_plot(1,:));
grid on
h=legend('$C_{switch}$', '$\overline{r}$', '$C_1$', '$C_2$', 'ed'); 
set(h,'Interpreter','latex')
plotParam

% 
% ax(2)=subplot(2,1,2);
% plot(t,J_optim,'b',t,J_normal, 'r');
% legend('Funzionale di Costo Reset', 'Funzionale di Costo Normal');   
% grid on
% plotParam   
% 
% linkaxes(ax,'x');
