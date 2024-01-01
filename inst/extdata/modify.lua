function raw_tex (t)
  return pandoc.RawBlock('tex', t)
end

--- Wrap code blocks in tcolorbox environments
function CodeBlock (cb)
  return {raw_tex'\\begin{snugshade}', cb, raw_tex '\\end{snugshade}'}
end
