function epBlock = defineBlocks(EEG)
% define trial blocks into IL, TR and PT
nep = length(EEG.epoch);
context = cell(nep, 1);
for i = 1:nep, context(i, 1) = EEG.epoch(i).eventcondType(1); end
[~, ~, tmp_block_id] = unique(context);
tmp_jump = diff(tmp_block_id);
tmp_b = 1;
b = 1;
for i_ep = 1:(nep - 1)
    if abs(tmp_jump(i_ep)) > 0
        block{b, 1} = tmp_b:i_ep;
        tmp_b = i_ep + 1;
        b = b + 1;
    end
end
block{b, 1} = tmp_b:nep;
% block ep into IL1, IL2~19, TR1, TR2~19, PT1, or PT2~57
epBlock = { block{1, 1}(:, 1), block{1, 1}(:, 2:end), ...
            block{2, 1}, [block{4:2:end, 1}], ...
            block{3, 1}(:, 1), [block{3, 1}(:, 2:end), block{5:2:end, 1}]};
end