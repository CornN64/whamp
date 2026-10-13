function ok = whamp_verify(refFile)
%WHAMP_VERIFY Re-dump all cases and diff against a stored baseline.
%   ok = WHAMP_VERIFY() compares against baseline.txt and prints the diff.
if nargin < 1, refFile = 'baseline.txt'; end
whamp_dump('current.txt');
a = fileread(refFile);
b = fileread('current.txt');
if strcmp(a, b)
    fprintf('IDENTICAL to %s\n', refFile);
    ok = true;
    return
end
ok = false;
la = strsplit(a, newline);
lb = strsplit(b, newline);
n = max(numel(la), numel(lb));
shown = 0;
for k = 1:n
    x = '';  if k <= numel(la), x = la{k}; end
    y = '';  if k <= numel(lb), y = lb{k}; end
    if ~strcmp(x, y)
        fprintf('L%d\n  base: %s\n  now : %s\n', k, x, y);
        shown = shown + 1;
        if shown >= 25
            fprintf('... (further differences suppressed)\n');
            return
        end
    end
end
