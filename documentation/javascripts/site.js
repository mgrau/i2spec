// Runs on every page: typeset the maths, and inline the figures so they take this page's theme colours.
function typeset(root) {
  if (window.renderMathInElement) {
    renderMathInElement(root, {
      delimiters: [{left: "\\[", right: "\\]", display: true}, {left: "\\(", right: "\\)", display: false}],
      throwOnError: false,
    });
  }
  for (const box of root.querySelectorAll("[data-svg]")) {
    if (box.dataset.done) continue;
    box.dataset.done = "1";
    fetch(box.dataset.svg).then(r => (r.ok ? r.text() : "")).then(t => { box.innerHTML = t; }).catch(() => {});
  }
}
if (window.document$) document$.subscribe(({body}) => typeset(body));
else addEventListener("load", () => typeset(document.body));
