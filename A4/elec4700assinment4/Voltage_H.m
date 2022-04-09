function [R3] = Voltage_H(vinlow,vinmax,limitSim)


nx = 200;
ny = 100;
limitSim = 100;

Curr = zeros(1,limitSim);
Necksizes = zeros(1,limitSim);
voltagesweep = linspace(vinlow, vinmax, limitSim);
for u = 1:limitSim

    Lb = floor(nx/3);
    Wb = floor(ny/3);
    Necksizes(u) = Wb;
    
    CuCond = 1;
    NoCond = 10e-2;
    
    %Conductivity map
    
    cMap = zeros(nx,ny);
    
    for i = 1:nx
        for j = 1: ny
            cMap(i,j) = CuCond;
        end
    end
    
    for i = 1:nx
        for j = 1:ny
            if (i>=1 && i<=Wb && j>Lb && j<=(2*Lb))
                cMap(i,j) = NoCond;
            end
            
            if (i<=ny && i>=(ny-Wb) && j>Lb && j<=(2*Lb))
                cMap(i,j) = NoCond;
            end
        end
    end
    
    G = sparse(nx*ny, nx*ny);
    F = zeros(1, nx*ny);
    
    for i = 1:nx
        for j = 1:ny
            n = j + (i - 1) * ny;
            
            if i == 1
                G(n,:) = 0;
                G(n,n) = 1;
                F(n) = voltagesweep(u);
                
            elseif i == nx
                G(n,:) = 0;
                G(n,n) = 1;
                
            elseif j == 1
                nxm = j + (i-2) * ny;
                nxp = j + (i) * ny;
                nyp = (j+1) + (i-1) * ny;
                
                rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
                rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
                ryp = (cMap(i,j) + cMap(i,j+1))/2.0;
                
                G(n,n) = -(rxm + rxp + ryp);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n, nyp) = ryp;
                
            elseif j == ny
                nxm = j + (i-2) * ny;
                nxp = j + (i) * ny;
                nym = (j-1) + (i-1) * ny;
                
                rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
                rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
                rym = (cMap(i,j) + cMap(i,j-1))/2.0;
                
                G(n,n) = -(rxm + rxp + rym);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n, nym) = rym;
                
            else
                nxm = j + (i-2) * ny;
                nxp = j + (i) * ny;
                nym = (j-1) + (i-1) * ny;
                nyp = (j+1) + (i-1) * ny;
                
                rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
                rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
                rym = (cMap(i,j) + cMap(i,j-1))/2.0;
                ryp = (cMap(i,j) + cMap(i,j+1))/2.0;
                
                G(n,n) = -(rxm + rxp + rym + ryp);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n, nym) = rym;
                G(n, nyp) = ryp;
            end
        end
    end
    V = G\F';
    
    for i = 1:nx
        for j = 1:ny
            n = j + (i-1) * ny;
            VG(i,j) = V(n);
        end
    end
    
    for i = 1:nx
        for j = 1:ny
            if i == 1
                Ex(i, j) = (VG(i + 1, j) - VG(i, j));
            elseif i == nx
                Ex(i, j) = (VG(i, j) - VG(i - 1, j));
            else
                Ex(i, j) = (VG(i + 1, j) - VG(i - 1, j)) * 0.5;
            end
            if j == 1
                Ey(i, j) = (VG(i, j + 1) - VG(i, j));
            elseif j == ny
                Ey(i, j) = (VG(i, j) - VG(i, j - 1));
            else
                Ey(i, j) = (VG(i, j + 1) - VG(i, j - 1)) * 0.5;
            end
        end
    end
    
    Ex = -Ex;
    Ey = -Ey;
    
    Jx = cMap .* Ex;
    Jy = cMap .* Ey;
    
    C0 = sum(Jx(1, :));
    Cnx = sum(Jx(nx, :));
    
    %Current
    Curr(u) = (C0 + Cnx) * 0.5;
    
    %Bottle Neck
    Necksizes(u) = ny - (2*Wb);
    
    %Sizes
    Wb = Wb+1;
end

%voltage = linspace(0.1, 30);
curr = linspace(0.1, 100);

currentvolt = polyfit(voltagesweep, Curr,1);

R3 = 1/currentvolt(1);
end