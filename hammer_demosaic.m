function [rgb_rec] = hammer_demosaic(res)

bayer_rec = single(zeros(512, 512, 22));
rgb_rec = single(zeros(512, 512, 3, 22));

for i=1:22
    bayer_rec(1:2:end, 1:2:end, :) = res(:, :, :, 1)*255;
    bayer_rec(1:2:end, 2:2:end, :) = res(:, :, :, 2)*255;
    bayer_rec(2:2:end, 1:2:end, :) = res(:, :, :, 3)*255;
    bayer_rec(2:2:end, 2:2:end, :) = res(:, :, :, 4)*255;
end

for i=1:22
    rgb_rec(:, :, :, i) = single(demosaic(uint16(reshape(bayer_rec(:, :, i), 512, 512)), 'rggb'))/255;
end

end