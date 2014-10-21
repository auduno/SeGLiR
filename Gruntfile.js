module.exports = function(grunt) {
  grunt.loadNpmTasks('grunt-contrib-concat');

  grunt.initConfig({
    concat: {
      dist: {
        src: [
            'source/head.js',
            'source/bernoulli.js',
            'source/tail.js'],
        dest: './SeGLiR.js' 
      }
    }
  });

  // Default task.
  grunt.registerTask('default', [
    'concat',
  ]);
};
