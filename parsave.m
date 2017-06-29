% Saves dataset, It is handy to save inside a parfor loop
function parsave(fname,dwpliPre,dwpliStg,dwpliTra,dwpliCIt,dwpliCIs,info,wpdTra,wpdStg,dwpliPreStd,dwpliPreMed,dwpliPreMad)
    save(fname, 'dwpliPre','dwpliStg','dwpliTra','dwpliCIt','dwpliCIs','info','wpdTra','wpdStg','dwpliPreStd','dwpliPreMed','dwpliPreMad')
end