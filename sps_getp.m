function v = sps_getp(s, name, dflt)
%SPS_GETP  从结构体 s 中安全取字段 name,不存在或为空时返回默认值 dflt
%
%   用于让各场模型只声明自己需要的参数,其余自动取默认值。

if isstruct(s) && isfield(s, name)
    v = s.(name);
    if isempty(v)
        v = dflt;
    end
else
    v = dflt;
end

end
