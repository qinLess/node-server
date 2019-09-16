const puppeteer = require('puppeteer')


class NewBrower {
    constructor(option) {
        this.option = option
        this.brower = null
    }

    async openBrower() {
        this.brower = await puppeteer.launch(this.option)
    }
}