module.exports = {
  env: { browser: true, es2021: true, node: true },
  extends: ['eslint:recommended', 'plugin:react/recommended'],
  parserOptions: { ecmaVersion: 'latest', sourceType: 'module' },
  settings: { react: { version: 'detect' } },
  rules: {
    'react/react-in-jsx-scope': 'off',
    // The project has no PropTypes dependency and performs no runtime prop
    // validation; component contracts are documented as JSDoc in
    // src/shared/analysisTypes.js instead.
    'react/prop-types': 'off',
    // Apostrophes and quotes in JSX text are intentional (mixed IT/EN copy).
    'react/no-unescaped-entities': 'off',
  },
}
