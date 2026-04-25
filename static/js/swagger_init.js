/**
 * swagger_init.js — boot SwaggerUIBundle for /docs.
 *
 * Lives in /static/ so the CSP (`script-src 'self' …`) admits it
 * without an `'unsafe-inline'` exception or per-request nonce/hash.
 *
 * The OpenAPI spec URL is fixed by the routes/openapi.py blueprint;
 * if the URL ever moves, expose it on the script tag via
 *   <script src="…/swagger_init.js" data-spec-url="…"></script>
 * and read currentScript.dataset.specUrl here.
 */
(function () {
  "use strict";
  function boot() {
    if (typeof SwaggerUIBundle === "undefined") return;
    window.ui = SwaggerUIBundle({
      url: "/api/v1/openapi.json",
      dom_id: "#swagger-ui",
      deepLinking: true,
      docExpansion: "list",
      tagsSorter: "alpha",
      operationsSorter: "alpha",
      tryItOutEnabled: true,
    });
  }
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", boot);
  } else {
    boot();
  }
})();
