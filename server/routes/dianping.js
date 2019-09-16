const router = require('koa-router')()
const dp = require('../layout/v2/dianping')


// 大众点评 图片token 加密
router.post('/dp_token', async (ctx, next) => {
    const { request } = ctx;
    let data = request.body;
    let token = dp.getToken(data.request_code, data.captcha_code, data.location)
    ctx.response.body = token
})

module.exports = router