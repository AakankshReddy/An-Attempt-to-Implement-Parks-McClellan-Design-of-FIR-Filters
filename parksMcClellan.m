w=linspace(0,pi,200000);
delta1=0.005;
delta2=0.005;
wp=0.4*pi;
ws=0.5*pi;
D=@(w)(w<=wp);
N=10;
np = round((N+2) * wp / (wp + pi - ws));   
ns = (N+2) - np;                            
omega=zeros(1,np+ns);
omega(1:np) = linspace(0, wp, np);
omega(np+1:end) = linspace(ws, pi, ns);
flag=1;
function [chk] = checkCONDITION(e, w, wp, ws, delta1, delta2)
chk = 0;
a = find(w <= wp, 1, 'last');   
b = find(w >= ws, 1, 'first'); 

peak1 = max(abs(e(1:a)));
peak2 = max(abs(e(b:end)));

if peak1 > delta1 || peak2 > delta2
    chk = 1;
end
end
%while 1
for j=1:1000
    signs = (-1).^(1:(N+2));
    cosMatrix = [cos(omega(:) * (0:N)), signs(:)];

    rhs = zeros(N+2,1);
    for i=1:np
        rhs(i) = 1;
    end
    for i=np+1:N+2
        rhs(i) = 0;
    end

    sol = pinv(cosMatrix) * rhs;
    g = sol(1:N+1);
    delta = sol(N+2);

    s=0;
    for (i=0:N)
        s = s + g(i + 1) * cos(i*w);
    end
    e=D(w)-s;

    if checkCONDITION(e, w, wp, ws, delta)
        [~, locsPos] = findpeaks(e);
        [~, locsNeg] = findpeaks(-e);
        locs = [locsPos, locsNeg];
        a = find(w <= wp, 1, 'last');
        b = find(w >= ws, 1, 'first');
        locs = unique([1, a, b, numel(w), locs]);
        minGap = round(numel(w) / (3*(N+2)));
        locs = sort(locs);
        keep = true(size(locs));
        lastKept = locs(1);
        for k = 2:numel(locs)
            if locs(k) - lastKept < minGap
                keep(k) = false;
            else
                lastKept = locs(k);
            end
        end
        locs = locs(keep);
        if numel(locs) > N+2
            forced = [1, a, b, numel(w)];
            rest = setdiff(locs, forced);
            [~, order] = sort(abs(e(rest)), 'descend');
            nRemaining = (N+2) - numel(forced);
            locs = sort([forced, rest(order(1:nRemaining))]);
        end
        if numel(locs) < N+2
            extra = (N+2) - numel(locs);
            fillLocs = round(linspace(1, numel(w), extra));
            locs = unique([locs, fillLocs]);
            locs = locs(1:(N+2));
        end
        omega = w(locs);
        np = sum(omega <= wp);
        ns = numel(omega) - np;
    else
        break
    end
end

figure;
plot(w, s, 'b', 'LineWidth', 1.2); 
