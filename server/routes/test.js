const router = require('koa-router')()

let code_msg = ''

router.get('/set_code', async (ctx, next) => {
    let { code } = ctx.request.query

    code_msg = code != '' && code ? code : ''

    ctx.response.body = code_msg
})

router.get('/get_code', async (ctx, next) => {
    // console.log('code_msg: ', code_msg)
    ctx.response.body = code_msg
    code_msg = ''
})

router.post('/test_post', async (ctx, next) => {
    ctx.response.body = '1111111111111'
})

module.exports = router