% ---------- Q1 CODE ----------
Ap = [0 1; -800 -1.6];
Bp = [0; 1];
Cp = [0 1]; % dy/dt is measured

%Q1.2. State Feedback Controller Design
CL = [Bp Ap*Bp]; %Controllability matrix
Pc = conv([1 3],[1 4]);
K = T2place(Ap, Bp, Pc);

%Q1.3. Observer Design
OL = [Cp; Cp*Ap]; %Observability matrix
Pob = conv([1 6],[1 7]);
Kob = T2place(Ap',Cp',Pob)';
