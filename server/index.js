const Koa = require('koa');
const app = new Koa();
const bodyparser = require('koa-bodyparser');
const static = require('koa-static')
const path = require('path')

const token = require('./routes/token');
const test = require('./routes/test');
const x_for_with = require('./routes/x_for_with');
const dp = require('./routes/dianping');

// Mozilla/5.0 (Linux; Android 6.0; Nexus 5 Build/MRA58N) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/72.0.3626.119 Mobile Safari/537.36

app.use(bodyparser());

const staticPath = './static'
app.use(static(
    path.join( __dirname, staticPath)
))

app.use(async (ctx, next) => {
  const start = new Date();
  await next();
  const ms = new Date() - start;
  console.log(`${ctx.method} ${ctx.url} - ${ms}ms`);
});

app.use(token.routes(), token.allowedMethods())
app.use(test.routes(), test.allowedMethods())
app.use(x_for_with.routes(), x_for_with.allowedMethods())
app.use(dp.routes(), dp.allowedMethods())

app.listen(80);
console.log('start server ! 127.0.0.1:80');

app.on('error', (err, ctx) => {
  console.error('server error', err, ctx)
});

module.exports = app;