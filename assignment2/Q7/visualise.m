%=====================================
% Plot 2D eigenvectors into a PDF file
%=====================================

try   
    U = load('solution.dat');
    n = size(U,2);
    m = sqrt(n);
    U = reshape(U', m, m, 2);
catch ME
    fprintf('Failed to load solution.dat: %s\n', ME.message)
    exit;
end

warning off MATLAB:hg:AutoSoftwareOpenGL

fprintf('Plotting the solution...\n');

x = (1:m)'/m;
y = repmat(x', m, 1);
x = repmat(x, 1, m);

subplot(2,1,1);
surf(x,y,U(:,:,1));
view(2);
colorbar;
xlabel('x')
ylabel('y')
title('u_1');

subplot(2,1,2);
surf(x,y,U(:,:,2));
view(2);
colorbar;
xlabel('x')
ylabel('y')
title('u_2');

p = get(gcf, 'Position');
set(gcf, 'Position', [0 0 p(3) p(4)*2]);

print(gcf,'-depsc2', 'solution.ps');
system('ps2pdf -dEPSCrop solution.ps');
delete('solution.ps');

fprintf('Success! solution.pdf was created\n');
exit;
