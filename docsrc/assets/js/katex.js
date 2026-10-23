// Renders the LaTeX formulas produced by pymdownx.arithmatex (also in the
// API docstrings), on first load and after every instant-navigation page change.
// KaTeX and its auto-render extension are loaded before this script (mkdocs.yml).
// Only the elements marked by arithmatex are rendered, so that a "$" or a
// parenthesis in normal text is never taken for a formula.
document$.subscribe(({ body }) => {
  if (typeof renderMathInElement !== "function") return
  body.querySelectorAll(".arithmatex").forEach((element) => {
    renderMathInElement(element, {
      delimiters: [
        { left: "\\(", right: "\\)", display: false },
        { left: "\\[", right: "\\]", display: true }
      ],
    })
  })
})
