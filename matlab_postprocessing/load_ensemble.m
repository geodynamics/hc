function [result] = load_ensemble(ensemble_file)
fh = fopen(ensemble_file,'r');
n = 0;
while ~feof(fh)
    line = fgetl(fh);
    if line(1) ~= '#'
        n = n + 1;
    end
end
fclose(fh);
chain = zeros(n,1);
residual = zeros(n,1);
nlayer = zeros(n,1);
var = zeros(n,1);
likeprob = zeros(n,1);
rad  = NaN*zeros(16,n);
visc = NaN*zeros(16,n);

i=1;
fh = fopen(ensemble_file,'r');
while ~feof(fh)
    line = fgetl(fh);
    if line(1) == '#'
        % do nothing
    else
        fields = sscanf(line,'%f,');
        if ~isempty(fields)
            chain(i) = fields(1);
            likeprob(i) = fields(4);
            residual(i) =fields(3);            
            var(i) = fields(5);
            nlayer(i) = fields(6);
            rad_start = 7;
            rad_end = rad_start + nlayer(i)-1;
            visc_start = rad_end+1;
            visc_end = visc_start + nlayer(i)-1;
            visc(1:nlayer(i),i) = fields(visc_start:visc_end);
            rad(1:nlayer(i),i) = fields(rad_start:rad_end);
            
            i = i + 1;
        else
            disp(['empty line at ' num2str(i)]);
        end
    end
    
end
result.residual = residual;
result.nlayer = nlayer;
result.var = var;
result.visc = visc;
result.rad = rad;
result.n = n;
result.chain = chain;