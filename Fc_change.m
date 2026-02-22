function sig = Fc_change(sig,freOffest)
%该版本中的freOffest为fc/fs
t_slot = 1:1:length(sig);
t_slot = reshape(t_slot,size(sig,1),size(sig,2));
sig = sig.*exp(1j*(2*pi*t_slot*freOffest));
end

