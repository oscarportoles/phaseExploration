% Saves dataset, It is handy to save inside a parfor loop
function parsave(fname, likehood, bumpMag, gammPara, eventprobs)
    save(fname, 'likehood', 'bumpMag', 'gammPara', 'eventprobs')
end